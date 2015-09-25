#!/usr/bin/env python
# Copyright (c) 2010 The George Washington University
# Authors: Chongzhi Zang, Weiqun Peng
#
# Disclaimer
# 
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
from GenomeData import *
import SeparateByChrom
import Utility

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def write(item, out):
	"""
	write one line into outfile. The file openning and closing is handled by outside. 
	item is a BED3 object
	"""
	#chrom, start, end, name, score, strand
	outline = item.chrom + "\t" + str(item.start) + "\t" + str(item.end) + "\n";
	out.write(outline);	


def union_islands_to_file(islandlist, f):
	islandlist.sort(key=operator.attrgetter('start'));
	current = islandlist[0]
	i = 1
	while i < len(islandlist):
		compare = islandlist[i]
		assert current.start <= compare.start
		if compare.start > current.end:
			write(current, f)
			current = compare
			i += 1
		else: 
			current.end = max(current.end, compare.end)
			i += 1
	write(current, f)


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--islandfile1", action="store", type="string", dest="islandfile1", metavar="<file>", help="file 1 with islands info to be unioned")
	parser.add_option("-b", "--islandfile2", action="store", type="string", dest="islandfile2", metavar="<file>", help="file 2 with islands info to be unioned; if no, type in any word")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8 or hg18", metavar="<str>")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	
	SeparateByChrom.separateByChrom(chroms, opt.islandfile1, '.island1')
	
	if Utility.fileExists(opt.islandfile2):
		SeparateByChrom.separateByChrom(chroms, opt.islandfile2, '.island2')
		
		for chrom in chroms: 
			f = open(chrom + '.output', 'w')
			bed_vals_1 = BED.BED(opt.species, chrom+'.island1', "BED3", 0)
			bed_vals_2 = BED.BED(opt.species, chrom+'.island2', "BED3", 0)
			if len(bed_vals_1[chrom]) > 0 or len(bed_vals_2[chrom]) > 0:
				islandlist = bed_vals_1[chrom] + bed_vals_2[chrom];
				union_islands_to_file(islandlist, f)
			f.close()
		SeparateByChrom.cleanup(chroms, '.island2')
	else:
		for chrom in chroms: 
			f = open(chrom + '.output', 'w')
			bed_vals_1 = BED.BED(opt.species, chrom+'.island1', "BED3", 0)
			if len(bed_vals_1[chrom]) > 0:
				islandlist = bed_vals_1[chrom]
				union_islands_to_file(islandlist, f)
			f.close()

	SeparateByChrom.combineAllGraphFiles(chroms, '.output', opt.outfile);
	
	SeparateByChrom.cleanup(chroms, '.output')
	SeparateByChrom.cleanup(chroms, '.island1')

if __name__ == "__main__":
	main(sys.argv)       