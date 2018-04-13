/*
 * @(#)TextViewer.java	1.11 00/05/24
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

import java.awt.*;
import java.io.*;
import java.beans.*;
import javax.activation.*;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JScrollPane;


/**
 * A very simple TextViewer Bean for the MIMEType "text/plain"
 *
 * @version	1.11, 00/05/24
 * @author	Christopher Cotton
 */

public class TextViewer extends JPanel implements CommandObject 
{

    private JTextArea text_area = null;
    private DataHandler dh = null;
    private String	verb = null;

    /**
     * Constructor
     */
    public TextViewer() {
	super(new GridLayout(1,1));

	// create the text area
	text_area = new JTextArea();
	text_area.setEditable(false);
	text_area.setLineWrap(true);

	// create a scroll pane for the JTextArea
	JScrollPane sp = new JScrollPane();
	sp.setPreferredSize(new Dimension(300, 300));
	sp.getViewport().add(text_area);
	
	add(sp);
    }


    public void setCommandContext(String verb, DataHandler dh)
	throws IOException {

	this.verb = verb;
	this.dh = dh;
	
	this.setInputStream( dh.getInputStream() );
    }


  /**
   * set the data stream, component to assume it is ready to
   * be read.
   */
  public void setInputStream(InputStream ins) {
      
      int bytes_read = 0;
      // check that we can actually read
      ByteArrayOutputStream baos = new ByteArrayOutputStream();
      byte data[] = new byte[1024];
      
      try {
	  while((bytes_read = ins.read(data)) >0)
		  baos.write(data, 0, bytes_read);
	  ins.close();
      } catch(Exception e) {
	  e.printStackTrace();
      }

      // convert the buffer into a string
      // place in the text area
      text_area.setText(baos.toString());

    }
}
