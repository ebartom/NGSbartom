#!/usr/bin/env python
# 
# Authors: Chongzhi Zang, Weiqun Peng
#
#This version 1) This version does not require either file's islands are non-overlapping. The output will be the same 
#as before (islands in file1 that are either overlapped or non-overlapped with islands in file2). i.e. if file1's islands 
#have overlap, the output could still have overlap. If file2's islands have overlap (which would be problematic in overlap 
#determination in the old version), this version will union them together before calculating the overlaps. 
#2) input or output files are not required to be exact BED format with certain number of columns. Together with the change 
#on BED.py, as long as the first 3 columns are chrom, start, end, this module will work, and keep all the information on the 
#rest columns in the output file. 

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
# Version 1.1  4/25/2011


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
import Utility
from GenomeData import *
import SeparateByChrom

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def are_islands_sorted(island_list):
	"""
	Check if sorted in ascending order.
	input is a list of BED with chrom, start, end and value.
	output: sorted =1 or 0
	"""
	sorted = 1;
	for index in range(0, len(island_list)-1):
		if island_list[index].start> island_list[index+1].start:
			sorted = 0;
	return sorted;


def union_islands(islandlist):
	"""
	The islandlist MUST be pre-sorted according to start!!!
	"""
	start_list =[]
	end_list = []
	current = islandlist[0]
	i = 1
	while i < len(islandlist):
		compare = islandlist[i]
		assert current.start <= compare.start
		if compare.start > current.end:
			start_list.append(current.start)
			end_list.append(current.end)
			current = compare
			i += 1
		else: 
			current.end = max(current.end, compare.end)
			i += 1
	start_list.append(current.start)
	end_list.append(current.end)
	return start_list, end_list


def region_overlap(start, end, start_list, end_list):
	"""
	Assuming non overlapping islands.
	"""
	assert (start <= end);
	start_position = bisect.bisect_right(end_list, start);
	end_position = bisect.bisect_left(start_list, end); 
	if (start_position < end_position): 
		return 1;
	else:
		return 0;


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--islandfile1", action="store", type="string", dest="islandfile1", metavar="<file>", help="file 1 with islands info to be compared")
	parser.add_option("-b", "--islandfile2", action="store", type="string", dest="islandfile2", metavar="<file>", help="file 2 with islands info to be compared")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8 or hg18", metavar="<str>")
	parser.add_option("-p", "--overlapin1", action="store", type="string", dest="overlapin1", metavar="<file>", help="file for islands in 1 overlapping with islands in 2")
	parser.add_option("-q", "--nonoverlapin1", action="store", type="string", dest="nonoverlapin1", help="file for islands in 1 not overlapping with islands in 2 ", metavar="<file>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	
	total_overlap_number_1 = 0
	total_islands_1 = 0
	
	SeparateByChrom.separateByChrom(chroms, opt.islandfile1, '.island1')
	SeparateByChrom.separateByChrom(chroms, opt.islandfile2, '.island2')
	
	for chrom in chroms: 
		f = open(chrom + '.1in2', 'w')
		g = open(chrom + '.1notin2', 'w')
		bed_vals_2 = BED.BED(opt.species, chrom+'.island2', "BED3", 0)
		if Utility.fileExists(chrom+'.island1') and len(bed_vals_2[chrom])>0:
			islandlist2 = bed_vals_2[chrom];
			if (are_islands_sorted(islandlist2) != 1):
				islandlist2.sort(key=operator.attrgetter('start'));
			(island2_start_list, island2_end_list) = union_islands(islandlist2)
			islands1 = open(chrom+'.island1', 'r')
			for line in islands1:
				if not re.match("#", line):
					total_islands_1 += 1
					line = line.strip()
					sline = line.split()
					start = int(sline[1])
					end = int(sline[2])
					if (region_overlap(start, end, island2_start_list, island2_end_list) == 1):
						f.write('\t'.join(sline) + '\n')
						total_overlap_number_1 += 1;
					else:
						g.write('\t'.join(sline) + '\n');
		elif Utility.fileExists(chrom+'.island1') and (len(bed_vals_2[chrom])==0):
			islands1 = open(chrom+'.island1', 'r')
			for line in islands1:
				if not re.match("#", line):
					total_islands_1 += 1
					line = line.strip()
					sline = line.split()
					g.write('\t'.join(sline) + '\n');
		f.close()
		g.close()
	
	print "total number of island in "+opt.islandfile1+":     ", total_islands_1;
	print "total number of island in "+opt.overlapin1+":     ", total_overlap_number_1;
	
	SeparateByChrom.combineAllGraphFiles(chroms, '.1in2', opt.overlapin1);
	SeparateByChrom.combineAllGraphFiles(chroms, '.1notin2', opt.nonoverlapin1);
	
	SeparateByChrom.cleanup(chroms, '.1in2')
	SeparateByChrom.cleanup(chroms, '.1notin2')
	SeparateByChrom.cleanup(chroms, '.island1')
	SeparateByChrom.cleanup(chroms, '.island2')

if __name__ == "__main__":
	main(sys.argv)       
