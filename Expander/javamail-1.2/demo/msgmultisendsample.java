/*
 * @(#)msgmultisendsample.java	1.12 00/05/24
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
import javax.mail.internet.*;
import javax.activation.*;

/**
 * msgmultisendsample creates a simple multipart/mixed message and sends it.
 * Both body parts are text/plain.
 * <p>
 * usage: <code>java msgmultisendsample <i>to from smtp true|false</i></code>
 * where <i>to</i> and <i>from</i> are the destination and
 * origin email addresses, respectively, and <i>smtp</i>
 * is the hostname of the machine that has smtp server
 * running.  The last parameter either turns on or turns off
 * debugging during sending.
 *
 * @author	Max Spivak
 */
public class msgmultisendsample {
    static String msgText1 = "This is a message body.\nHere's line two.";
    static String msgText2 = "This is the text in the message attachment.";

    public static void main(String[] args) {
	if (args.length != 4) {
	    System.out.println("usage: java msgmultisend <to> <from> <smtp> true|false");
	    return;
	}

	String to = args[0];
	String from = args[1];
	String host = args[2];
	boolean debug = Boolean.valueOf(args[3]).booleanValue();

	// create some properties and get the default Session
	Properties props = new Properties();
	props.put("mail.smtp.host", host);

	Session session = Session.getDefaultInstance(props, null);
	session.setDebug(debug);
	
	try {
	    // create a message
	    MimeMessage msg = new MimeMessage(session);
	    msg.setFrom(new InternetAddress(from));
	    InternetAddress[] address = {new InternetAddress(to)};
	    msg.setRecipients(Message.RecipientType.TO, address);
	    msg.setSubject("JavaMail APIs Multipart Test");
	    msg.setSentDate(new Date());

	    // create and fill the first message part
	    MimeBodyPart mbp1 = new MimeBodyPart();
	    mbp1.setText(msgText1);

	    // create and fill the second message part
	    MimeBodyPart mbp2 = new MimeBodyPart();
	    // Use setText(text, charset), to show it off !
	    mbp2.setText(msgText2, "us-ascii");

	    // create the Multipart and its parts to it
	    Multipart mp = new MimeMultipart();
	    mp.addBodyPart(mbp1);
	    mp.addBodyPart(mbp2);

	    // add the Multipart to the message
	    msg.setContent(mp);
	    
	    // send the message
	    Transport.send(msg);
	} catch (MessagingException mex) {
	    mex.printStackTrace();
	    Exception ex = null;
	    if ((ex = mex.getNextException()) != null) {
		ex.printStackTrace();
	    }
	}
    }
}
