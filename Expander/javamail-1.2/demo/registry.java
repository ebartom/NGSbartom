/*
 * @(#)registry.java	1.8 00/05/24
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
import java.util.*;

import javax.mail.*;
import javax.mail.internet.*;

/**
 * This class demonstrates how to query the registry for available
 * Providers, set default providers, etc. See section 5.2 in the
 * JavaMail 1.0 Spec for details on how to use the registry.
 * 
 * See the comments inline for what's happening.
 *
 * @author Max Spivak
 */

public class registry {
    // let's remember a few providers
    static Provider _aProvider, _bProvider, _sunSMTP, _sunIMAP; 

    public static void main(String[] args) {
	Properties props = new Properties();

	// set smtp and imap to be our default 
	// transport and store protocols, respectively
	props.put("mail.transport.protocol", "smtp");
	props.put("mail.store.protocol", "imap");

	// 
	props.put("mail.smtp.class", "com.sun.mail.smtp.SMTPTransport");
	props.put("mail.imap.class", "com.sun.mail.imap.IMAPStore");
	
	Session session = Session.getInstance(props, null);
	//session.setDebug(true);

	// Retrieve all configured providers from the Session
	System.out.println("\n------ getProviders()----------");
	Provider[] providers = session.getProviders();
	for (int i = 0; i < providers.length; i++) {
	    System.out.println("** " + providers[i]);

	    // let's remember some providers so that we can use them later
	    // (I'm explicitly showing multiple ways of testing Providers)
	    // BTW, no Provider "ACME Corp" will be found in the default
	    // setup b/c its not in any javamail.providers resource files
	    String s = null;
	    if (((s = providers[i].getVendor()) != null) &&
		s.startsWith("ACME Corp")) {
		_aProvider = providers[i];
	    }

	    // this Provider won't be found by default either
	    if (providers[i].getClassName().endsWith("application.smtp"))
		_bProvider = providers[i];

	    // this Provider will be found since com.sun.mail.imap.IMAPStore
	    // is configured in javamail.default.providers
	    if (providers[i].getClassName().equals("com.sun.mail.imap.IMAPStore")){
		_sunIMAP = providers[i];
	    }

	    // this Provider will be found as well since there is a
	    // Sun Microsystems SMTP transport configured by 
	    // default in javamail.default.providers
	    if (((s = providers[i].getVendor()) != null) &&
		s.startsWith("Sun Microsystems") && 
		providers[i].getType() == Provider.Type.TRANSPORT &&
		providers[i].getProtocol().equalsIgnoreCase("smtp")) {
		_sunSMTP = providers[i];
	    }
	}
	
	System.out.println("\n------ initial protocol defaults -------");
	try {
	    System.out.println("imap: " + session.getProvider("imap"));
	    System.out.println("smtp: " + session.getProvider("smtp"));
	    // the NNTP provider will fail since we don't have one configured
	    System.out.println("nntp: " + session.getProvider("nntp"));
	} catch (NoSuchProviderException mex) { 
	    System.out.println(">> This exception is OK since there is no NNTP Provider configured by default");
	    mex.printStackTrace();
	}

	System.out.println("\n------ set new protocol defaults ------");
	// set some new defaults
	try {
	    // since _aProvider isn't configured, this will fail
	    session.setProvider(_aProvider);  // will fail
	} catch (NoSuchProviderException mex) { 
	    System.out.println(">> Exception expected: _aProvider is null");
	    mex.printStackTrace(); 
	}
	try {
	    // _sunIMAP provider should've configured correctly; should work
	    session.setProvider(_sunIMAP);
	} catch (NoSuchProviderException mex) { mex.printStackTrace(); }
	try {
	    System.out.println("imap: " + session.getProvider("imap"));
	    System.out.println("smtp: " + session.getProvider("smtp"));
	} catch (NoSuchProviderException mex) { mex.printStackTrace(); }


	System.out.println("\n\n----- get some stores ---------");
	// multiple ways to retrieve stores. these will print out the
	// string "imap:" since its the URLName for the store
	try {
	    System.out.println("getStore(): " + session.getStore());
	    System.out.println("getStore(Provider): " + 
			       session.getStore(_sunIMAP));
	} catch (NoSuchProviderException mex) {
	    mex.printStackTrace();
	}

	try {
	    System.out.println("getStore(imap): " + session.getStore("imap"));
	    // pop3 will fail since it doesn't exist
	    System.out.println("getStore(pop3): " + session.getStore("pop3"));
	} catch (NoSuchProviderException mex) { 
	    System.out.println(">> Exception expected: no pop3 provider");
	    mex.printStackTrace(); 
	}


	System.out.println("\n\n----- now for transports/addresses ---------");
	// retrieve transports; these will print out "smtp:" (like stores did)
	try {
	    System.out.println("getTransport(): " + session.getTransport());
	    System.out.println("getTransport(Provider): " + 
			       session.getTransport(_sunSMTP));
	    System.out.println("getTransport(smtp): " + 
			       session.getTransport("smtp"));
	    System.out.println("getTransport(Address): " + 
			       session.getTransport(new InternetAddress("mspivak@apilon")));
	    // News will fail since there's no news provider configured
	    System.out.println("getTransport(News): " + 
			       session.getTransport(new NewsAddress("rec.humor")));
	} catch (MessagingException mex) { 
	    System.out.println(">> Exception expected: no news provider configured");
	    mex.printStackTrace(); 
	}
    }
}
