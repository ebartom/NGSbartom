/*
 * @(#)StoreTreeNode.java	1.8 00/05/24
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

import javax.swing.tree.DefaultMutableTreeNode;
import javax.mail.*;

/**
 * Node which represents a Store in the javax.mail apis. 
 *
 * @version 1.8, 00/05/24
 * @author Christopher Cotton
 */
public class StoreTreeNode extends DefaultMutableTreeNode {
    
    protected Store	store = null;
    protected Folder	folder = null;
    protected String	display = null;

    /**
     * creates a tree node that points to the particular Store.
     *
     * @param what	the store for this node
     */
    public StoreTreeNode(Store what) {
	super(what);
	store = what;
    }

    
    /**
     * a Store is never a leaf node.  It can always contain stuff
     */
    public boolean isLeaf() {
	return false;
    }
   

    /**
     * return the number of children for this store node. The first
     * time this method is called we load up all of the folders
     * under the store's defaultFolder
     */

    public int getChildCount() {
	if (folder == null) {
	    loadChildren();
	}
	return super.getChildCount();
    }
    
    protected void loadChildren() {
	try {
	    // connect to the Store if we need to
	    if (!store.isConnected()) {
		store.connect();
	    }

	    // get the default folder, and list the
	    // subscribed folders on it
	    folder = store.getDefaultFolder();
	    // Folder[] sub = folder.listSubscribed();
	    Folder[] sub = folder.list();

	    // add a FolderTreeNode for each Folder
	    int num = sub.length;
	    for(int i = 0; i < num; i++) {
		FolderTreeNode node = new FolderTreeNode(sub[i]);
		// we used insert here, since add() would make
		// another recursive call to getChildCount();
		insert(node, i);
	    }
	    
	} catch (MessagingException me) {
	    me.printStackTrace();
	}
    }

    /**
     * We override toString() so we can display the store URLName
     * without the password.
     */

    public String toString() {
	if (display == null) {
	    URLName url = store.getURLName();
	    if (url == null) {
		display = store.toString();
	    } else {
		// don't show the password
		URLName too = new URLName( url.getProtocol(), url.getHost(), url.getPort(),
					   url.getFile(), url.getUsername(), null);
		display = too.toString();
	    }
	}
	
	return display;
    }
    
    
}

