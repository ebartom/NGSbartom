/*
 * @(#)msgsendsample.java	1.17 00/06/14
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
 * msgsendsample creates a very simple text/plain message and sends it.
 * <p>
 * usage: <code>java msgsendsample <i>to from smtphost true|false</i></code>
 * where <i>to</i> and <i>from</i> are the destination and
 * origin email addresses, respectively, and <i>smtphost</i>
 * is the hostname of the machine that has the smtp server
 * running. The last parameter either turns on or turns off
 * debugging during sending.
 *
 * @author Max Spivak
 */
public class msgsendsample {
    static String msgText = "This is a message body.\nHere's the second line.";

    public static void main(String[] args) {
	if (args.length != 4) {
	    usage();
	    System.exit(1);
	}

	System.out.println();
	
	String to = args[0];
	String from = args[1];
	String host = args[2];
	boolean debug = Boolean.valueOf(args[3]).booleanValue();

	// create some properties and get the default Session
	Properties props = new Properties();
	props.put("mail.smtp.host", host);
	if (debug) props.put("mail.debug", args[3]);

	Session session = Session.getDefaultInstance(props, null);
	session.setDebug(debug);
	
	try {
	    // create a message
	    Message msg = new MimeMessage(session);
	    msg.setFrom(new InternetAddress(from));
	    InternetAddress[] address = {new InternetAddress(args[0])};
	    msg.setRecipients(Message.RecipientType.TO, address);
	    msg.setSubject("JavaMail APIs Test");
	    msg.setSentDate(new Date());
	    // If the desired charset is known, you can use
	    // setText(text, charset)
	    msg.setText(msgText);
	    
	    Transport.send(msg);
	} catch (MessagingException mex) {
	    System.out.println("\n--Exception handling in msgsendsample.java");

	    mex.printStackTrace();
	    System.out.println();
	    Exception ex = mex;
	    do {
		if (ex instanceof SendFailedException) {
		    SendFailedException sfex = (SendFailedException)ex;
		    Address[] invalid = sfex.getInvalidAddresses();
		    if (invalid != null) {
			System.out.println("    ** Invalid Addresses");
			if (invalid != null) {
			    for (int i = 0; i < invalid.length; i++) 
				System.out.println("         " + invalid[i]);
			}
		    }
		    Address[] validUnsent = sfex.getValidUnsentAddresses();
		    if (validUnsent != null) {
			System.out.println("    ** ValidUnsent Addresses");
			if (validUnsent != null) {
			    for (int i = 0; i < validUnsent.length; i++) 
				System.out.println("         "+validUnsent[i]);
			}
		    }
		    Address[] validSent = sfex.getValidSentAddresses();
		    if (validSent != null) {
			System.out.println("    ** ValidSent Addresses");
			if (validSent != null) {
			    for (int i = 0; i < validSent.length; i++) 
				System.out.println("         "+validSent[i]);
			}
		    }
		}
		System.out.println();
		if (ex instanceof MessagingException)
		    ex = ((MessagingException)ex).getNextException();
		else
		    ex = null;
	    } while (ex != null);
	}
    }

    private static void usage() {
	System.out.println("usage: java msgsendsample <to> <from> <smtp> true|false");
    }
}
