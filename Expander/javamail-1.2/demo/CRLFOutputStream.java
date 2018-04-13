/*
 * @(#)CRLFOutputStream.java	1.2 00/05/24
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

import java.io.*;

/**
 * Convert lines into the canonical MIME format, that is,
 * terminate lines with CRLF. <p>
 *
 * This stream can be used with the Part.writeTo and Message.writeTo
 * methods to generate the canonical MIME format of the data for the
 * purpose of (e.g.) sending it via SMTP or computing a digital
 * signature.
 */
public class CRLFOutputStream extends FilterOutputStream {
    protected int lastb = -1;
    protected static byte[] newline;
    static {
	newline = new byte[2];
	newline[0] = (byte)'\r';
	newline[1] = (byte)'\n';
    }

    public CRLFOutputStream(OutputStream os) {
	super(os);
    }

    public void write(int b) throws IOException {
	if (b == '\r') {
	    out.write(newline);
	} else if (b == '\n') {
	    if (lastb != '\r')
		out.write(newline);
	} else {
	    out.write(b);
	}
	lastb = b;
    }

    public void write(byte b[]) throws IOException {
	write(b, 0, b.length);
    }

    public void write(byte b[], int off, int len) throws IOException {
	int start = off;

	len += off;
	for (int i = start; i < len ; i++) {
	    if (b[i] == '\r') {
		out.write(b, start, i - start);
		out.write(newline);
		start = i + 1;
	    } else if (b[i] == '\n') {
		if (lastb != '\r') {
		    out.write(b, start, i - start);
		    out.write(newline);
		}
		start = i + 1;
	    }
	    lastb = b[i];
	}
	if ((len - start) > 0)
	    out.write(b, start, len - start);
    }
}
