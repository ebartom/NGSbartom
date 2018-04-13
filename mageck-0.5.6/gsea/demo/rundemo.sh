#!/bin/bash

# g++ -I./ -O3 -o calculategsea calculategsea.cpp

# time ./calculategsea -r ~/Dropbox/work/crispr/linearmodel/permu/next9.abl_veh_plasmid.txt -g ../c2.cp.kegg.v4.0.symbols.gmt -p 10000

# time ./calculategsea -r ~/Dropbox/work/crispr/linearmodel/permu/next9.abl_veh_plasmid.txt -g ../kegg.ribosome.gmt -p 10000

time  ../../bin/mageckGSEA -r demo1.txt -g kegg.ribosome.gmt -p 100000 
