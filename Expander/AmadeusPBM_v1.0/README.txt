========================================
AMADEUS v1.2 (Dec 2008)
ALLEGRO v1.0 (Dec 2008)
AMADEUSPBM V1.0 (Jul 2011)
========================================

AMADEUS - A Motif Algorithm for Detecting Enrichment in MUltiple Species.
ALLEGRO - A Log-Likelihood based Engine for Gene expression Regulatory motifs Over-representation discovery.
AMADEUSPBM - A Motif finding algorithm for extracting binding site motifs from PBM data

AMADEUS and ALLEGRO are developed at Ron Shamir's Computational Genomics Group, Tel Aviv University.
Copyright: Tel Aviv University, 2007-2011.

The Amadeus/Allegro/AmadeusPBM software is freely available for academic use under the conditions
of a Tel Aviv University licence agreement.


Publications:
-------------
Allegro:
"Allegro: Analyzing expression and sequence in concert to discover regulatory programs",
Y. Halperin, C. Linhart, I. Ulitsky and Ron Shamir,
Nucleic Acids Research, to appear.

Amadeus:
"Transcription factor and microRNA motif discovery: The Amadeus platform and a compendium of metazoan target sets",
C. Linhart, Y. Halperin and R. Shamir,
Genome Research, vol. 18:7, pp. 1180-1189, July 2008.


Web-site:
---------
Allegro: http://acgt.cs.tau.ac.il/allegro (also for Amadeus and AmadeusPBM)
Amadeus: http://acgt.cs.tau.ac.il/amadeus
The above web-sites contain information about the software (documentation,
algorithm, results, etc.), as well as downloadable sequence files (promoters
and 3'UTRs from various organisms).
Please visit the web-sites for important updates.


Installation and execution
--------------------------
NOTE: Amadeus/Allegro/AmadeusPBM requires Java version 1.5.
The AmadeusPBM_v1.0.zip file should be extracted into the directory in which you wish
to place Allegro (the same software contains Amadeus and AmadeusPBM).
It will automatically generate a directory called "AmadeusPBM_v1.0".
Run the file "run.bat" (double-click on it in the Windows file explorer) in order
to launch Allegro. This should open the Allegro window. The left panel in this
window is used to enter the input data and parameters; the right panel shows the
output results.
Set the "Data type" parameter according to the type of analysis you wish to perform:
a. "Target set" - search for motifs that are over-represented in the promoters (or 3' UTRs)
of a given list of genes. This type of analysis is done using our Amadeus algorithm.
b. "Expression" - search for motifs that have a distinct expression profile, given one
or more gene expression datasets. This type of analysis is done using our Allegro algorithm.
c. "PBM" - search for motifs in PBM data using the AmadeusPBM algorithm.

Supplied are three sample datasets, one for each type of analysis:

a. Target-set analysis:
In order to run Amadeus on the sample data supplied with the installation, first
download the human promoters zip file from http://acgt.cs.tau.ac.il/allegro/download.html.
Next, click on the right button at the bottom of the left panel ("Load parameters from file"),
and choose the file "params/params_HCC.txt" - this will load the parameters for the
sample human cell-cycle data. This target set consists of 344 genes that peak
in the G2 and G2/M phases [Whitfield et al., 2002]. Next, click on the left
button ("Run") to start the execution. Amadeus should now run on the sample data,
and print the textual output to the "Output" tab in the right panel. The final
results are shown in the "Results" tab. Using the supplied parameters, the top
scoring motifs Amadeus detects in the G2+G2/M data are CHR (TTTGAA or TTCAAA)
and NF-Y (CCAAT box).

b. Expression analysis:
In order to run Allegro on the sample data supplied with the installation, first
download the mouse promoters zip file from http://acgt.cs.tau.ac.il/allegro/download.html.
Next, click on the right button at the bottom of the left panel ("Load parameters from file"),
and choose the file "params/params_TLRs.txt" - this will load the parameters for the
sample TLRs (Toll-lie receptors) expression data. This gene expression dataset was
downloaded from the Innate-Immunity System-Biology project (http://www.systemsbiology-immunity.org).
It contains expression profiles of RAW264.7 monocyte macrophage-like cell line after
exposure to several pathogen-mimetic agents. Each of the agents used is recognized by a
different subset of TLRs. Exposure to lipopolysaccharide (LPS) was sampled at several
time points.
Next, click on the left button ("Run") to start the execution. Allegro should now run
on the sample data, and print the textual output to the "Output" tab in the right panel.
The final results are shown in the "Results" tab. Using the supplied parameters, the top
scoring motifs Allegro detects are ISRE/IRF (AA..GAAACT or AGTTTC..TT)
and NF-kB (GG[AG]..T[CT]CC).

c. PBM analysis:
In order to run AmadeusPBM on the sample data supplied with the installation,
click on the right button at the bottom of the left panel ("Load parameters from file"),
and choose the file "params/params_PBM.txt" - this will load the parameters for the
sample PBM data. These data were downloaded from Uniprobe database
(http://the_brain.bwh.harvard.edu/uniprobe/index.php), which contains PBM files of many proteins.
The file format that AmadeusPBM accepts is _deBruijn.txt.
Next, click on the left button ("Run") to start the execution. AmadeusPBM should now run
on the sample data, and print the textual output to the "Output" tab in the right panel.
The final results are shown in the "Results" tab. Using the supplied parameters,
the top scoring motifs AmadeusPBM detects are E2F, both primary and secondary motifs.

Files and directories
---------------------
AmadeusPBM_v1.0 directory contains the following:
- "run.bat" file - used to execute the program.
- "run_1.3G_mem.bat" file - used to execute the program with a larger memory,
e.g., when analyzing multiple datasets from several organisms, or when running
in "Large" mode. Note that this execution requires at least 2GB RAM.
- "AmadeusPBM_v1.0.jar" and "third_party.jar" files - Java archives with the
class files for running Amadeus/Allegro/AmadeusPBM.
- "licenses" directory - contains licensing information.
- "data/TFs" directory - used to store TRANSFAC matrix files with
position weight matrixes of known transcription factors. The motifs
found de-novo by Amadeus/Allegro/AmadeusPBM are compared to the TRANSFAC matrixes.
The Allegro installation contains the public TRANSFAC matrix table (v7.0,
30 Sep 2005).
- "data/miRNA" directory - used to store matrix files with micro-RNA
binding seeds. The installation contains 8-bases seeds downloaded from
miRBase (release 16.0, Sep 2010). See also data/miRNA/README.txt.
- "data/sequences" directory - used to store promoter and 3' UTR sequences.
Sequences for various organisms, as well as the corresponding BG sets,
can be downloaded from http://acgt.cs.tau.ac.il/allegro/download.html.
- "targetSets" directory - contains a sample target-set file with human
G2+G2/M cell-cycle genes.
- "expression" directory - contains a sample gene expression file with
the mouse TLRs dataset.
- "params" directory - contains parameters files for running Amadeus/Allegro/AmadeusPBM
on the sample datasets.
- "PBM" directory - contains a pair of smaple PBM files of TF E2F2.


Minimum hardware and system requirements
----------------------------------------
- O/S: Amadeus/Allegro/AmadeusPBM runs under Windows (2000/XP) and Linux.
Allegro was tested on Windows XP professional (SP2), and on Linux Debian.
Allegro should also run on most Unix flavors, and other OS's that support Java 1.5.
- CPU: Pentium 4, or equivalent
- RAM memory size: 1GB
Analyzing multiple datasets, or running in mode "Large", requires 2GB RAM (use the run_1.3G_mem.bat file);
If this doesn't suffice, please edit the bat file, and set a larger memory limit (the -Xmx parameter),
according to your system.
- Disk space: The Allegro program requires 20 MB of hard-disk space,
and additional space for the sequence files (download from http://acgt.cs.tau.ac.il/amadeus).
- Java VM: Java version 1.5


License
-------
The Amadeus/Allegro/AmadeusPBM software is freely available for academic use.
It is also available for non-academic use under appropriate licensing.
Please contact Ron Shamir (rshamir@post.tau.ac.il),
Chaim Linhart (chaiml@post.tau.ac.il) or Yaron Orenstien (yaronore@post.tau.ac.il)
for further information.

Amadeus/Allegro/AmadeusPBM uses several third-party Java libraries that are distributed under
various license agreements, as detailed in the file licenses/licenses.txt.


Main changes in AmadeusPBM v1.0:
--------------------------------
- Implemented the AmadeusPBM algorithm for extracing TFBS from PBM data.

Main changes in Allegro v1.0:
-----------------------------
- Implemented the Allegro algorithm for gene expression analysis.
- Changed default values of various parameters.

Main changes in Amadeus v1.2:
-----------------------------
- New step in the "greedy PWM" phase for enumerating small modifications to the PWM.
- In 3' UTR analysis, the discovered motifs are compared only to the rev-comp strand of the miRNA seeds.

Main changes in Amadeus v1.1:
-----------------------------
Added an option to run Amadeus on user-supplied sequences.


Contact info
------------
In case you encounter problems, or if you have any questions on Amadeus/Allegro/AmadeusPBM,
please contact us:
Chaim Linhart (chaiml@post.tau.ac.il)
Yaron Orenstien (yaronore@post.tau.ac.il)
Yonit Halperin (yonithalp@gmail.com)
Igor Ulitsky (ulitskyi@post.tau.ac.il)
