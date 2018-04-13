/*
 * @(#)SimpleClient.java	1.23 00/05/24
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
/*
 * @(#)SimpleClient.java	1.23, 00/05/24
 *
 * Copyright (c) 1997-1998 by Sun Microsystems, Inc.
 * All Rights Reserved.
 */

import java.util.*;
import java.io.*;
import javax.mail.*;
import javax.mail.internet.*;
import javax.activation.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.table.*;
import javax.swing.tree.*;
import javax.swing.event.*;


/**
 * Demo app that shows a very simple Mail Client
 *
 * @version 1.23, 00/05/24
 * @author Christopher Cotton
 * @author Bill Shannon
 */

public class SimpleClient {

    static Vector url = new Vector();
    static FolderViewer fv;
    static MessageViewer mv;

    public static void main(String argv[]) {
	boolean usage = false;

	for (int optind = 0; optind < argv.length; optind++) {
	    if (argv[optind].equals("-L")) {
		url.addElement(argv[++optind]);
	    } else if (argv[optind].startsWith("-")) {
		usage = true;
		break;
	    } else {
		usage = true;
		break;
	    }
	}

	if (usage || url.size() == 0) {
	    System.out.println("Usage: SimpleClient -L url");
	    System.out.println("   where url is protocol://username:password@hostname/");
	    System.exit(1);
	}

        try {
	    // Set up our Mailcap entries.  This will allow the JAF
	    // to locate our viewers.
	    File capfile = new File("simple.mailcap");
	    if (!capfile.isFile()) {
		System.out.println(
		    "Cannot locate the \"simple.mailcap\" file.");
		System.exit(1);
	    }
	    
	    CommandMap.setDefaultCommandMap( new MailcapCommandMap(
		new FileInputStream(capfile)));
		
 	    JFrame frame = new JFrame("Simple JavaMail Client");
 	    frame.addWindowListener(new WindowAdapter() {
 		public void windowClosing(WindowEvent e) {System.exit(0);}});
	    //frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
	    // Get a Store object
	    SimpleAuthenticator auth = new SimpleAuthenticator(frame);
	    Session session = 
		Session.getDefaultInstance(System.getProperties(), auth);
	    //session.setDebug(true);

 	    DefaultMutableTreeNode root = new DefaultMutableTreeNode("Root");

	    // create a node for each store we have
	    for (Enumeration e = url.elements() ; e.hasMoreElements() ;) {
		String urlstring = (String) e.nextElement();
		URLName urln = new URLName(urlstring);
		Store store = session.getStore(urln);
		
		StoreTreeNode storenode = new StoreTreeNode(store);
		root.add(storenode);
	    }	    

	    DefaultTreeModel treeModel = new DefaultTreeModel(root);
	    JTree tree = new JTree(treeModel);
	    tree.addTreeSelectionListener(new TreePress());

	    /* Put the Tree in a scroller. */
	    JScrollPane        sp = new JScrollPane();
	    sp.setPreferredSize(new Dimension(250, 300));
	    sp.getViewport().add(tree);

	    /* Create a double buffered JPanel */
	    JPanel sv = new JPanel(new BorderLayout());
	    sv.add("Center", sp);

	    fv = new FolderViewer(null);

	    JSplitPane jsp = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
				sv, fv);
	    jsp.setOneTouchExpandable(true);
	    mv = new MessageViewer();
	    JSplitPane jsp2 = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
				jsp, mv);
	    jsp2.setOneTouchExpandable(true);

	    frame.getContentPane().add(jsp2);
	    frame.pack();
	    frame.show();

	} catch (Exception ex) {
	    System.out.println("SimpletClient caught exception");
	    ex.printStackTrace();
	    System.exit(1);
	}
    }

}

class TreePress implements TreeSelectionListener {
    
    public void valueChanged(TreeSelectionEvent e) {
	TreePath path = e.getNewLeadSelectionPath();
	if (path != null) {
	    Object o = path.getLastPathComponent();
	    if (o instanceof FolderTreeNode) {
		FolderTreeNode node = (FolderTreeNode)o;
		Folder folder = node.getFolder();

		try {
		    if ((folder.getType() & Folder.HOLDS_MESSAGES) != 0) {
			SimpleClient.fv.setFolder(folder);
		    }
		} catch (MessagingException me) { }
	    }
	}
    }
}
