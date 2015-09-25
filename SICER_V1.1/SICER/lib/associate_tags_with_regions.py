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
import bisect

import BED
import UCSC
import GenomeData;
import Utility
import SeparateByChrom
import get_total_tag_counts

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();

def tag_position(sline, fragment_size):
	shift = int(round(fragment_size/2))
	if plus.match(sline[5]):
		return atoi(sline[1]) + shift
	elif minus.match(sline[5]):
		return atoi(sline[2]) - 1 - shift

def countTagsInWindow(start, end, tag_starts):
	# Require that the tag_starts are sorted!
	assert( start<=end )
	start_ind = bisect.bisect_left(tag_starts, start);
	end_ind = bisect.bisect_right(tag_starts, end);
	tags = end_ind - start_ind;
	return tags;


def find_readcount_on_regions(tag_position_list, region_start_list, region_end_list):
	'''
	The regions could overlap !!
	returns a list, with total tag number on this region, order as the region lists
	'''
	region_readcount_list = []
	assert len(region_start_list) == len(region_end_list)
	if (Utility.is_list_sorted(tag_position_list)==0):
		tag_position_list.sort();
	for i in range(0, len(region_start_list)):
		tags = countTagsInWindow(region_start_list[i], region_end_list[i], tag_position_list);
		region_readcount_list.append(tags);
	return region_readcount_list

def find_readcount_on_islands(island_start_list, island_end_list, tag):
	"""
	Make sure the islands are sorted.
	Islands are non-overlapping!
	Returns the index of the island on which the tag lands, or -1.
	"""
	
	index = bisect.bisect_right(island_start_list, tag);
	if index - bisect.bisect_left(island_end_list, tag) == 1:
		return index-1;
	else:
		return -1;
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--rawreadfile", action="store", type="string", dest="readfile", metavar="<file>", help="raw read file in bed format")
	parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", metavar="<int>", help="average size of a fragment after CHIP experiment")
	parser.add_option("-b", "--islandfile", action="store", type="string", dest="islandfile", metavar="<file>", help="island file in BED format")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="island read count file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	
	
	islands = BED.BED(opt.species, opt.islandfile, "BED3", 0);
	if Utility.fileExists(opt.readfile):
		SeparateByChrom.separateByChrom(chroms, opt.readfile, '.bed1');
	else:
		print opt.readfile, " not found";
		sys.exit(1)
	
	total = 0; 
	library_size = get_total_tag_counts.get_total_tag_counts(opt.readfile);
	
	scaling_factor = 1000000; 
	out = open(opt.out_file, 'w');
	for chrom in chroms:
		if chrom in islands.keys():
			island_list = islands[chrom];
			island_readcount_list=[0]*len(island_list);
			
			if Utility.is_bed_sorted(island_list) == 0:
				island_list.sort(key=operator.attrgetter('start'));
				
			island_start_list = []
			island_end_list = []
			for item in island_list:
				island_start_list.append(item.start)
				island_end_list.append(item.end)

			read_file = chrom + ".bed1";
			f = open(read_file,'r')
			for line in f:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					position = tag_position(sline, opt.fragment_size)
					index = find_readcount_on_islands(island_start_list, island_end_list, position);
					if index >= 0:
						island_readcount_list[index] += 1;
						total += 1;
			f.close();
							
			
			for index in xrange(len(island_list)):
				item = island_list[index];
				normalized_read_count = island_readcount_list[index]/float(library_size) * scaling_factor;
				outline = item.chrom + "\t" + str(item.start) + "\t" + str(item.end) + "\t" + str(island_readcount_list[index]) +  "\t" + str(normalized_read_count) + "\n";	
				out.write(outline);		
							
	SeparateByChrom.cleanup(chroms, '.bed1');
	out.close();
	print "Total number of reads on islands are: ", total; 


if __name__ == "__main__":
	main(sys.argv)
