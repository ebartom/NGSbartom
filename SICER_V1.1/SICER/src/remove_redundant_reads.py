#!/usr/bin/env python
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

import GenomeData
import SeparateByChrom
import Utility


def remove_redundant_1chrom_single_strand_sorted(infile, outfile, cutoff):
	'''infile can only contain reads from one chromosome and only one kind of strands (+/-). file must be pre-sorted by column2, then by column3.'''
	f = open(infile,'r')
	o = open(outfile, 'w')
	current_start = 0
	current_end = 0
	current_count = 1
	total = 0
	retained = 0
	for line in f:
		if not re.match("#", line):
			total += 1
			line = line.strip()
			sline = line.split()
			start = atoi(sline[1])
			end = atoi(sline[2])
			if start != current_start:
				o.write('\t'.join(sline)+'\n')
				retained += 1
				current_start = start
				current_end = end
				current_count = 1
			elif end != current_end:
				o.write('\t'.join(sline)+'\n')
				retained += 1
				current_start = start
				current_end = end
				current_count = 1
			else:
				current_count += 1
				assert current_start == start
				assert current_end == end
				if current_count <= cutoff:
					o.write('\t'.join(sline)+'\n')
					retained += 1
	f.close()
	o.close()
	return (total, retained)


def strand_broken_remove(chrom, cutoff):
	'''infile can only contain reads from one chromosome'''
	infile = chrom + ".bed1";
	outfile = chrom + ".bed2"
	try:
		if os.system('grep [[:space:]]+ %s | sort -g -k 2,3 > plus.bed1' % (infile)):
			raise
	except: sys.stderr.write("+ reads do not exist in " + str(infile) + "\n");
	(p_total, p_retained) = remove_redundant_1chrom_single_strand_sorted('plus.bed1', 'plus_removed.bed1', cutoff)
	
	try:
		if os.system('grep [[:space:]]- %s | sort -g -k 2,3 > minus.bed1' % (infile)):
			raise
	except: sys.stderr.write("- reads do not exist in " + str(infile) + "\n");
	(m_total, m_retained) = remove_redundant_1chrom_single_strand_sorted('minus.bed1', 'minus_removed.bed1', cutoff)
	
	print chrom, "\tPlus reads:",p_total, "\tRetained plus reads:", p_retained,     ";\tMinus reads:", m_total, "\tRetained minus reads:", m_retained;
	
	os.system('cat plus_removed.bed1 minus_removed.bed1 > %s' % (outfile))
	os.system('rm plus*.bed1')
	os.system('rm minus*.bed1')


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="species under consideration", metavar="<str>")
	parser.add_option("-b", "--raw_bed_file", action="store", type="string",
                      dest="bed_file", help="raw bed file", metavar="<file>")
	parser.add_option("-t", "--threshold", action="store", type="int",
                      dest="threshold", help="threshold for copy number", metavar="<int>")	      
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	SeparateByChrom.separateByChrom(chroms, opt.bed_file, '.bed1')
	
	for chrom in chroms:
		if (Utility.fileExists(chrom + ".bed1")):
			strand_broken_remove(chrom, opt.threshold)
	
	SeparateByChrom.combineAllGraphFiles(chroms, '.bed2', opt.out_file)
	SeparateByChrom.cleanup(chroms, '.bed1')
	SeparateByChrom.cleanup(chroms, '.bed2')


if __name__ == "__main__":
	main(sys.argv)
