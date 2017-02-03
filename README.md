![alt text](https://github.com/ebartom/NGSbartom/blob/master/CetoDiagram2.png "Ceto")

Ceto is a modular system for analyzing next generation sequence (NGS) data.  Within Ceto are many R scripts and perl scripts (the Ceto Toolbox) which are shared in this directory.

buildPipelineScripts.pl is essentially a decision tree that walks through the tools in the toolbox along with open source tools, and creates a series of shell scripts for different types of NGS analysis.  Currently, ChIP-seq and RNA-seq pipelines are implemented and well tested.  These are based on collaborations between the Bartom Lab, the Shilatifard lab, the Ntziachristos Lab, the Pulmonology Division, and other scientists at Northwestern University.

To read more about the philosophy behind Ceto, please see https://github.com/ebartom/NGSbartom/blob/master/000_Introduction.txt .

buildPipelineScripts.pl should be stable, while buildPipelineScripts.draft.pl will contain more experimental analyses.

If you find a bug in Ceto, or just have an improvement to suggest, please email ebartom@northwestern.edu with Ceto in the subject line.
