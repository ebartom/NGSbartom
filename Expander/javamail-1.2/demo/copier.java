/*
 * @(#)copier.java	1.7 00/05/24
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

/**
 *
 * @version	1.7, 00/05/24
 * @author	Christopher Cotton
 */

import javax.mail.*;

/**
 * copier will copy a specified number of messages from one folder
 * to another folder. it demonstrates how to use the JavaMail APIs
 * to copy messages.<p>
 *
 * usage for copier: copier <i>protocol</i> <i>host</i> <i>user</i> 
 * <i>password</i> <i>src folder</i> <i>dest folder</i> <i>start msg #</i> <i>end msg #</i><p>
 *
 */

public class copier {

  public static void main(String argv[]) {
      boolean debug = false;// change to get more errors
      
      if (argv.length != 5) {
	  System.out.println( "usage: copier <urlname> <src folder>" +
			      "<dest folder> <start msg #> <end msg #>");
	  return;
      }

      try {
	  URLName url = new URLName(argv[0]);
	  String src = argv[1];	// source folder
	  String dest = argv[2];	// dest folder
	  int start = Integer.parseInt(argv[3]);  // copy from message #
	  int end = Integer.parseInt(argv[4]);	// to message #

	  // Get the default Session object

	  Session session = Session.getDefaultInstance(
	      System.getProperties(), null);
	  // session.setDebug(debug);

	  // Get a Store object that implements 
	  // the protocol.
	  Store store = session.getStore(url);
	  store.connect();
	  System.out.println("Connected...");

	  // Open Source Folder
	  Folder folder = store.getFolder(src);
	  folder.open(Folder.READ_WRITE);
	  System.out.println("Opened source...");	  

	  if (folder.getMessageCount() == 0) {
		System.out.println("Source folder has no messages ..");
		folder.close(false);
		store.close();
	  }

	  // Open destination folder, create if needed 
	  Folder dfolder = store.getFolder(dest);
	  if (!dfolder.exists()) // create
	      dfolder.create(Folder.HOLDS_MESSAGES);
	  System.out.println("Opened dest...");	  

	  Message[] msgs = folder.getMessages(start, end);
	  System.out.println("Got messages...");	  

	  // Copy messages into destination, 
	  folder.copyMessages(msgs, dfolder);
	  System.out.println("Copied messages...");	  

	  // Close the folders and store
	  folder.close(false);
	  dfolder.close(false);
	  store.close();
	  System.out.println("Closed folders and store...");
	  
      } catch (Exception e) {
	  e.printStackTrace();
      }

      System.exit(0);
  }

}


