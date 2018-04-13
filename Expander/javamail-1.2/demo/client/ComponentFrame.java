/*
 * @(#)ComponentFrame.java	1.7 00/05/24
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
import java.awt.event.*;
import javax.swing.JFrame;
import javax.swing.WindowConstants;


/**
 * this Frame provides a utility class for displaying a single
 * Component in a Frame.
 *
 * @version	1.7, 00/05/24
 * @author	Christopher Cotton
 */

public class ComponentFrame extends JFrame {
    
    /**
     * creates the frame
     * @param what	the component to display
     */
    public ComponentFrame(Component what) {
	this(what, "Component Frame");
    }

    /**
     * creates the frame with the given name
     * @param what	the component to display
     * @param name	the name of the Frame
     */
    public ComponentFrame(Component what, String name) {
	super(name);

	// make sure that we close and dispose ourselves when needed
	setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);

	// default size of the frame
	setSize(700,600);

	// we want to display just the component in the entire frame
	if (what != null) {
	    getContentPane().add("Center", what);
	}
    }
}
