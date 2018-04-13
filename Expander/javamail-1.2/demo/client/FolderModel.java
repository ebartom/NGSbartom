/*
 * @(#)FolderModel.java	1.12 00/11/20
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
import java.util.Date;
import javax.swing.table.AbstractTableModel; 

/**
 * Maps the messages in a Folder to the Swing's Table Model
 *
 * @version	1.12, 00/11/20
 * @author	Christopher Cotton
 * @author	Bill Shannon
 */

public class FolderModel extends AbstractTableModel {
    
    Folder	folder;
    Message[]	messages;

    String[]	columnNames = { "Date", "From", "Subject"}; 
    Class[]	columnTypes = { String.class, String.class, String.class }; 

    public void setFolder(Folder what) throws MessagingException {
	if (what != null) {

	    // opened if needed
	    if (!what.isOpen()) {
		what.open(Folder.READ_WRITE);
	    }
    
	    // get the messages
	    messages = what.getMessages();
	    cached = new String[messages.length][];
	} else {
	    messages = null;
	    cached = null;
	}
	// close previous folder and switch to new folder
	if (folder != null)
	    folder.close(true);
	folder = what;
	fireTableDataChanged();
    }
    
    public Message getMessage(int which) {
	return messages[which];
    }

    //---------------------
    // Implementation of the TableModel methods
    //---------------------

    public String getColumnName(int column) {
	return columnNames[column];
    }
    
    public Class getColumnClass(int column) {
	return columnTypes[column];
    }
    

    public int getColumnCount() {
        return columnNames.length; 
    }

    public int getRowCount() {
	if (messages == null)
	    return 0;
	
	return messages.length;
    }
 
    public Object getValueAt(int aRow, int aColumn) {
	switch(aColumn) {
	case 0:	// date
	case 1: // From		String[] what = getCachedData(aRow);
	case 2: // Subject
	    String[] what = getCachedData(aRow);
	    if (what != null) {
		return what[aColumn];
	    } else {
		return "";
	    }
	    
	default:
	    return "";
	}
    }

    protected static String[][]	cached;
    
    protected String[] getCachedData(int row) {
	if (cached[row] == null) {
	    try{
		Message m = messages[row];
	    
		String[] theData = new String[4];
	    
		// Date
		Date date = m.getSentDate();
		if (date == null) {
		    theData[0] = "Unknown";
		} else {
		    theData[0] = date.toString();
		}
	    
		// From
		Address[] adds = m.getFrom();
		if (adds != null && adds.length != 0) {
		    theData[1] = adds[0].toString();	    
		} else {
		    theData[1] = "";
		}
		
		// Subject
		String subject = m.getSubject();
		if (subject != null) {
		    theData[2] = subject;
		} else {
		    theData[2] = "(No Subject)";
		}

		cached[row] = theData;
	    }
	    catch (MessagingException e) {
		e.printStackTrace();
	    }
	}
	
	return cached[row];
    }
}
