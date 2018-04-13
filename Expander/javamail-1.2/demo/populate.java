/*
 * @(#)populate.java	1.6 00/05/24
 *
 * Copyright 1997-2000 Sun Microsystems, Inc. All Rights Reserved.
 *
 * Sun grants you ("Licensee") a non-exclusive, royalty free, license to use,
 * modify and redistribute this software in source and binary code form,
 * provided that i) this copyright notice and license appear on all copies of
 * the software; and ii) Licensee does not utilize the software in a manner
 * which is disparaging to Sun.
 *
 * This software is provided "AS IS," without a warranty of any kind. ALL
 * EXPRESS OR IMPLIED CONDITIONS, REPRESENTATIONS AND WARRANTIES, INCLUDING ANY
 * IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE OR
 * NON-INFRINGEMENT, ARE HEREBY EXCLUDED. SUN AND ITS LICENSORS SHALL NOT BE
 * LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE AS A RESULT OF USING, MODIFYING
 * OR DISTRIBUTING THE SOFTWARE OR ITS DERIVATIVES. IN NO EVENT WILL SUN OR ITS
 * LICENSORS BE LIABLE FOR ANY LOST REVENUE, PROFIT OR DATA, OR FOR DIRECT,
 * INDIRECT, SPECIAL, CONSEQUENTIAL, INCIDENTAL OR PUNITIVE DAMAGES, HOWEVER
 * CAUSED AND REGARDLESS OF THE THEORY OF LIABILITY, ARISING OUT OF THE USE OF
 * OR INABILITY TO USE SOFTWARE, EVEN IF SUN HAS BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGES.
 *
 * This software is not designed or intended for use in on-line control of
 * aircraft, air traffic, aircraft navigation or aircraft communications; or in
 * the design, construction, operation or maintenance of any nuclear
 * facility. Licensee represents and warrants that it will not use or
 * redistribute the Software for such purposes.
 */

import javax.mail.*;
import javax.mail.internet.*;

/*
 * Copy folder hierarchies between different Stores. This is a useful 
 * utility to populate new (and possibly empty) mail stores. Specify
 * both the source and destination folders as URLs.
 *	
 * @author John Mani
 */

public class populate {

    static boolean force = false;
    static boolean skipSCCS = false;
    static boolean clear = false;

    public static void main(String argv[]) {
    	String srcURL = null;
    	String dstURL = null;
	boolean debug = false;

	int optind;

	for (optind = 0; optind < argv.length; optind++) {
	    if (argv[optind].equals("-s")) {
		srcURL = argv[++optind];
	    } else if (argv[optind].equals("-d")) {
		dstURL = argv[++optind];
	    } else if (argv[optind].equals("-D")) {
		debug = true;
	    } else if (argv[optind].equals("-f")) {
		force = true;
	    } else if (argv[optind].equals("-S")) {
		skipSCCS = true;
	    } else if (argv[optind].equals("-c")) {
		clear = true;
	    } else if (argv[optind].equals("--")) {
		optind++;
		break;
	    } else if (argv[optind].startsWith("-")) {
		printUsage();
		System.exit(1);
	    } else {
		break;
	    }
	}

	try {

	    if (srcURL == null || dstURL == null) {
		printUsage();
		System.exit(1);
	    }

	    Session session = Session.getDefaultInstance(
				System.getProperties(), null);
	    session.setDebug(debug);

	    // Get source folder
	    Folder srcFolder = session.getFolder(new URLName(srcURL));
	    if (!srcFolder.exists()) {
		System.out.println("source folder does not exist");
		srcFolder.getStore().close();
		System.exit(1);
	    }

	    // Set up destination folder
	    URLName dstURLName = new URLName(dstURL);
	    Folder dstFolder;
	    // Check if the destination URL has a folder specified. If
	    // not, we use the source folder name
	    if (dstURLName.getFile() == null) {
		Store s = session.getStore(dstURLName);
		s.connect();
		dstFolder = s.getFolder(srcFolder.getName());
	    } else
		dstFolder = session.getFolder(new URLName(dstURL));

	    if (clear && dstFolder.exists()) {
		if (!dstFolder.delete(true)) {
		    System.out.println("couldn't delete " +
						dstFolder.getFullName());
		    return;
		}
	    }
	    copy(srcFolder, dstFolder);

	    // Close the respective stores.
	    srcFolder.getStore().close();
	    dstFolder.getStore().close();

	} catch (MessagingException mex) {
	    System.out.println(mex.getMessage());
	    mex.printStackTrace();
	}
    }

    private static void copy(Folder src, Folder dst)
		throws MessagingException {
	System.out.println("Populating " + dst.getFullName());

	if (!dst.exists()) {
	    // Create it.
	    if (!dst.create(src.getType())) {
		System.out.println("couldn't create " + dst.getFullName());
		return;
	    }

	    // Copy over any messges from src to dst
	    if ((src.getType() & Folder.HOLDS_MESSAGES) != 0) {
		src.open(Folder.READ_ONLY);
		src.copyMessages(src.getMessages(), dst);
		src.close(false);
	    }
	} else  {
	    System.out.println(dst.getFullName() + " already exists");
	    // Copy over any messges from src to dst
	    if (force && (src.getType() & Folder.HOLDS_MESSAGES) != 0) {
		src.open(Folder.READ_ONLY);
		src.copyMessages(src.getMessages(), dst);
		src.close(false);
	    }
	}

	// Copy over subfolders
	if ((src.getType() & Folder.HOLDS_FOLDERS) != 0) {
		Folder[] sf = src.list();
		for (int  i = 0; i < sf.length; i++) {
		    // skip SCCS directories?
		    if (skipSCCS && sf[i].getName().equals("SCCS"))
			continue;
		    copy(sf[i], dst.getFolder(sf[i].getName()));
		}
    	}
    }

    private static void printUsage() {
	System.out.println("populate [-D] [-f] [-S] [-c] " +
			   "-s source_url -d dest_url");
	System.out.println("URLs are of the form: " +
		  	   "protocol://username:password@hostname/foldername");
	System.out.println("The destination URL does not need a foldername," +
		  	   " in which case, the source foldername is used");
    }
}
