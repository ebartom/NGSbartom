/*
 * @(#)ByteArrayDataSource.java	1.2 00/05/24
 *
 * Copyright 1998-2000 Sun Microsystems, Inc. All Rights Reserved.
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

import java.io.*;
import javax.activation.*;

/**
 * A simple DataSource for demonstration purposes.
 * This class implements a DataSource from:
 * 	an InputStream
 *	a byte array
 * 	a String
 *
 * @author John Mani
 * @author Bill Shannon
 * @author Max Spivak
 */
public class ByteArrayDataSource implements DataSource {
    private byte[] data;	// data
    private String type;	// content-type

    /* Create a DataSource from an input stream */
    public ByteArrayDataSource(InputStream is, String type) {
        this.type = type;
        try { 
            ByteArrayOutputStream os = new ByteArrayOutputStream();
	    int ch;

	    while ((ch = is.read()) != -1)
                // XXX - must be made more efficient by
	        // doing buffered reads, rather than one byte reads
	        os.write(ch);
	    data = os.toByteArray();

        } catch (IOException ioex) { }
    }

    /* Create a DataSource from a byte array */
    public ByteArrayDataSource(byte[] data, String type) {
        this.data = data;
	this.type = type;
    }

    /* Create a DataSource from a String */
    public ByteArrayDataSource(String data, String type) {
	try {
	    // Assumption that the string contains only ASCII
	    // characters!  Otherwise just pass a charset into this
	    // constructor and use it in getBytes()
	    this.data = data.getBytes("iso-8859-1");
	} catch (UnsupportedEncodingException uex) { }
	this.type = type;
    }

    /**
     * Return an InputStream for the data.
     * Note - a new stream must be returned each time.
     */
    public InputStream getInputStream() throws IOException {
	if (data == null)
	    throw new IOException("no data");
	return new ByteArrayInputStream(data);
    }

    public OutputStream getOutputStream() throws IOException {
	throw new IOException("cannot do this");
    }

    public String getContentType() {
        return type;
    }

    public String getName() {
        return "dummy";
    }
}
