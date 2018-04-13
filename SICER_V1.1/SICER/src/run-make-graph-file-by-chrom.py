#!/usr/bin/env python
# Authors: Dustin E Schones, Chongzhi Zang, Weiqun Peng and Keji Zhao
# Disclaimer
# 
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

## get BED module
import BED
import GenomeData
import make_graph_file
import SeparateByChrom

def makeGraphFile(chroms, chrom_lengths, window, fragment_size):
	for chrom in chroms:
		if chrom in chrom_lengths.keys():
			chrom_length = chrom_lengths[chrom];
		else:
			 print "Can not find the length of ", chrom;
			 exit(1);
		bed_file = chrom + ".bed";	
		graph_file = chrom + ".graph";
		make_graph_file.make_graph_file(bed_file, chrom, chrom_length, window, fragment_size, graph_file)		
 
def main(argv):
    """
    Note the window_size and the fragment_size are both input as strings, as they are used in
    a shell script in makeGraphFile.
    """
    parser = OptionParser()
    parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="mm8,hg18,dm2,etc", metavar="<str>")
    parser.add_option("-b", "--bed_file", action="store", type="string",
                      dest="bedfile", help="bed file to make graph file of",
                      metavar="<file>")
    parser.add_option("-w", "--window_size", action="store", type="int",
                      dest="window_size", help="window size", metavar="<int>")
    parser.add_option("-i", "--fragment_size", action="store", type="int",
                      dest="fragment_size",
                      help="size of fragments after CHIP experiment",
                      metavar="<int>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output bed summary file name",
                      metavar="<file>")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 10:
        parser.print_help()
        sys.exit(1)

    if opt.species in GenomeData.species_chroms.keys():
        chroms = GenomeData.species_chroms[opt.species];
	chrom_lengths = GenomeData.species_chrom_lengths[opt.species];
        SeparateByChrom.separateByChrom(chroms, opt.bedfile, ".bed");
	makeGraphFile(chroms, chrom_lengths, opt.window_size, opt.fragment_size);
        final_output_file = opt.outfile;
        final_output_file = SeparateByChrom.combineAllGraphFiles(chroms, ".graph", final_output_file);
        SeparateByChrom.cleanup(chroms, ".bed");
	SeparateByChrom.cleanup(chroms, ".graph");
    else:
        print opt.species + " is not in the species list ";
	
    

if __name__ == "__main__":
	main(sys.argv)


        
