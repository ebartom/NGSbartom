===========================================
 EXPANDER v6 
===========================================

EXPANDER - EXPression ANalyzer and DisplayER is a gene-expression analysis and visualization software.
EXPANDER is developed at Ron Shamir's Computational Genomics Group, Tel Aviv University.
Copyright: Tel Aviv University, 2003.
The Expander web-site is at: http://acgt.cs.tau.ac.il/expander/


Execution:
----------
Double click on the file "Expander.bat" (OR on "Expander_2GB.bat" if RAM size>=2GB OR "Expander_4GB.bat" if RAM size>=4GB ) that is located under the "Expander" directory.

When running on Linux/Unix OS:
- Make sure that you have "rwx" permissions for the Expander directory and for the directory in which your data is located.
- Make sure that you have "rx" permissions for all "*.exe" files that are under your Expander directory.


Files and directories:
----------------------
Expander directory contains the following:
 - "Expander.bat" file - used to execute the program (by double clicking it).
 - "Expander_2GB.bat" file - used to execute the program  (by double clicking it) in computers with at list 2GB RAM (for better performance).
 - "Expander_4GB.bat" file - used to execute the program  (by double clicking it) in computers with at list 4GB RAM (for better performance).
 - "help.html" file - a user help manual.
 - "ReleseNotes.txt" file - contains information regarding the changes made since the previous version.
 - "updates.html" file - contains direct access to the Expander download page.
 - "BSD.txt","MIT_caryoscope.txt", "LGPL.txt" and "copyrights.txt" files - contain copyright and licensing information.
 - "organisms" directory - used to store all organism-specific data (annotation files and TF fingerprint files,
    downloaded from the Expander download page should be extracted into this directory).
 - "sample_input_files" directory - contains sample data files that can be loaded into Expander.
 - Additional files & directories for internal usage of Expander.


Minimum hardware and system requirements:
-----------------------------------------------------------------------
- O/S: Expander runs under Windows (2000/XP) and Linux.
  Expander was tested on Windows XP professional (SP2), and on Linux Debian.
  Expander should also run on most Unix flavors.
- CPU: Pentium 4, or equivalent
- RAM memory size >= 1GB.
* If RAM size >= 2GB, Expander can be un using "Expander_2GB.bat" for better performance.
* If RAM size >= 4B, Expander can be un using "Expander_4GB.bat" for better performance.
- Disk space: The Expander program requires 16 MB of hard-disk space.
  Functional/Promoter analysis require additional 15-350MB per organism, depending
  on the organism.
- Java VM: Java version 5 or later.

Additional softwares required:
----------------------------------------------
R: The .cel file preprocessing and the newly added SAM filter utilities require the pre-installation of one of the recent versions of R, a free software environment for statistical computing and graphics. R can be installed from: http://cran.r-project.org/. 
After installing R, please do the following to install the Bioconductor “affy” package and the "samr" package:

1.      Run R.
2.      In the R frame\window type the text: source("http://bioconductor.org/biocLite.R")
3.      Press ‘Enter’.
4.      In the R frame\window type the text:  biocLite("affy")
5.      Press ‘Enter’.

To install the 'samr' package:
6.      In the R frame\window type the text:  install.packages("samr")
7.      Press ‘Enter’.

To install the 'ISA' package:
8.      In the R frame\window type the text: source("http://www.bioconductor.org/biocLite.R")
9.      Press ‘Enter’.
10.  In the R frame\window type the text biocLite("eisa")
11.  Press ‘Enter’.

You may install only one of the packages, depending on what you wish to use (to install only "samr", follow instructions number 1, 6 and 7).


Contact info:
-------------
In case you encounter problems, or if you have any questions on Expander,
please use the following email address:
expander@cs.tau.ac.il
