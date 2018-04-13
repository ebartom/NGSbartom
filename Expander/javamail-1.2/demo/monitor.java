/*
 * @(#)monitor.java	1.6 00/05/24
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

import java.util.*;
import java.io.*;
import javax.mail.*;
import javax.mail.event.*;
import javax.activation.*;

/* Monitors given mailbox for new mail */

public class monitor {

    public static void main(String argv[])
    {
	if (argv.length != 5) {
  	   System.out.println("Usage: monitor <host> <user> <password> <mbox> <freq>");
	   System.exit(1);
	}
	System.out.println("\nTesting monitor\n");

        try {
	    Properties props = System.getProperties();

	    // Get a Session object
	    Session session = Session.getDefaultInstance(props, null);
	    // session.setDebug(true);

	    // Get a Store object
	    Store store = session.getStore("imap");

	   // Connect
	    store.connect(argv[0], argv[1], argv[2]);

	   // Open a Folder
	    Folder folder = store.getFolder(argv[3]);
	    if (folder == null || !folder.exists()) {
		System.out.println("Invalid folder");
		System.exit(1);
	    }

	    folder.open(Folder.READ_WRITE);

	    // Add messageCountListener to listen for new messages
	    folder.addMessageCountListener(new MessageCountAdapter() {
		public void messagesAdded(MessageCountEvent ev) {
		    Message[] msgs = ev.getMessages();
		    System.out.println("Got " + msgs.length + " new messages");

		    // Just dump out the new messages
		    for (int i = 0; i < msgs.length; i++) {
			try {
			    DataHandler dh = msgs[i].getDataHandler();
			    InputStream is = dh.getInputStream();
			    int c;
			    while ((c = is.read()) != -1)
			 	System.out.write(c); 
			} catch (IOException ioex) { 
			    ioex.printStackTrace();	
			} catch (MessagingException mex) {
			    mex.printStackTrace();
			}
		    }
		}
	    });
			
	   // Check mail once in "freq" MILLIseconds

	    int freq = Integer.parseInt(argv[4]);

	    for (; ;) {
		Thread.sleep(freq); // sleep for freq milliseconds

		// This is to force the IMAP server to send us
		// EXISTS notifications. 
		folder.getMessageCount();
	    }

	} catch (Exception ex) {
	    ex.printStackTrace();
	}
    }
}
