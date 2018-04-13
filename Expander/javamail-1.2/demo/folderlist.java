/*
 * @(#)folderlist.java	1.27 00/08/30
 *
 * Copyright 1996-2000 Sun Microsystems, Inc. All Rights Reserved.
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

import java.util.Properties;
import javax.mail.*;

/**
 * Demo app that exercises the Message interfaces.
 * List information about folders.
 *
 * @author John Mani
 * @author Bill Shannon
 */

public class folderlist {
    static String protocol = null;
    static String host = null;
    static String user = null;
    static String password = null;
    static String url = null;
    static String root = null;
    static String pattern = "%";
    static boolean recursive = false;
    static boolean verbose = false;
    static boolean debug = false;

    public static void main(String argv[]) throws Exception {
	int optind;
	for (optind = 0; optind < argv.length; optind++) {
	    if (argv[optind].equals("-T")) {
		protocol = argv[++optind];
	    } else if (argv[optind].equals("-H")) {
		host = argv[++optind];
	    } else if (argv[optind].equals("-U")) {
		user = argv[++optind];
	    } else if (argv[optind].equals("-P")) {
		password = argv[++optind];
	    } else if (argv[optind].equals("-L")) {
		url = argv[++optind];
	    } else if (argv[optind].equals("-R")) {
		root = argv[++optind];
	    } else if (argv[optind].equals("-r")) {
		recursive = true;
	    } else if (argv[optind].equals("-v")) {
		verbose = true;
	    } else if (argv[optind].equals("-D")) {
		debug = true;
	    } else if (argv[optind].equals("--")) {
		optind++;
		break;
	    } else if (argv[optind].startsWith("-")) {
		System.out.println(
"Usage: folderlist [-T protocol] [-H host] [-U user] [-P password] [-L url]");
		System.out.println(
"\t[-R root] [-r] [-v] [-D] [pattern]");
		System.exit(1);
	    } else {
		break;
	    }
	}
	if (optind < argv.length)
	    pattern = argv[optind];

	// Get a Properties object
	Properties props = System.getProperties();

	// Get a Session object
	Session session = Session.getDefaultInstance(props, null);
	session.setDebug(debug);

	// Get a Store object
	Store store = null;
	Folder rf = null;
	if (url != null) {
	    URLName urln = new URLName(url);
	    store = session.getStore(urln);
	    store.connect();
	} else {
	    if (protocol != null)
		store = session.getStore(protocol);
	    else
		store = session.getStore();

	    // Connect
	    if (host != null || user != null || password != null)
		store.connect(host, user, password);
	    else
		store.connect();
	}

	// List namespace
	if (root != null)
	    rf = store.getFolder(root);
	else
	    rf = store.getDefaultFolder();

	dumpFolder(rf, false, "");
	if ((rf.getType() & Folder.HOLDS_FOLDERS) != 0) {
	    Folder[] f = rf.list(pattern);
	    for (int i = 0; i < f.length; i++)
		dumpFolder(f[i], recursive, "    ");
	}

	store.close();
    }

    static void dumpFolder(Folder folder, boolean recurse, String tab)
					throws Exception {
    	System.out.println(tab + "Name:      " + folder.getName());
	System.out.println(tab + "Full Name: " + folder.getFullName());
	System.out.println(tab + "URL:       " + folder.getURLName());

	if (verbose) {
	    if (!folder.isSubscribed())
		System.out.println(tab + "Not Subscribed");

	    if ((folder.getType() & Folder.HOLDS_MESSAGES) != 0) {
		if (folder.hasNewMessages())
		    System.out.println(tab + "Has New Messages");
		System.out.println(tab + "Total Messages:  " +
						folder.getMessageCount());
		System.out.println(tab + "New Messages:    " +
						folder.getNewMessageCount());
		System.out.println(tab + "Unread Messages: " +
						folder.getUnreadMessageCount());
	    }
	    if ((folder.getType() & Folder.HOLDS_FOLDERS) != 0)
		System.out.println(tab + "Is Directory");
	}

	System.out.println();

	if ((folder.getType() & Folder.HOLDS_FOLDERS) != 0) {
	    if (recurse) {
		Folder[] f = folder.list();
		for (int i = 0; i < f.length; i++)
		    dumpFolder(f[i], recurse, tab + "    ");
	    }
	}
    }
}
