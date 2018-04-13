#!/usr/bin/env python
# 
# Authors: Weiqun Peng
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
import random
import shutil
import get_total_tag_counts


def slice(desired_number_tags, raw_bed_file, out_file_name):
	"""
	Read a raw bed file and take thedesired number of lines
	
	"""
	current_total = get_total_tag_counts.get_total_tag_counts(raw_bed_file);
	if current_total <= desired_number_tags:
		# copy the file.
		shutil.copy(raw_bed_file, out_file_name); 
		print "existing number of tags is ", current_total,  "<= desired, no need to sample";
	else:
		count =0.0;
		infile = open(raw_bed_file,'r');
		outfile = open(out_file_name, 'w');
		for line in infile:
			if count >= desired_number_tags: break;	
			outfile.write(line);
			count +=1;
		outfile.close();
		infile.close();

		
def random_sample (desired_number_tags, raw_bed_file, out_file_name):
	"""
	Read a raw bed file and take the desired number of lines
	
	If reproduceable result is needed, one needs to set seed of the random function outside.
	
	"""
	current_total = int(get_total_tag_counts.get_total_tag_counts(raw_bed_file));
	assert  (current_total >= desired_number_tags);
	sample_list = random.sample(xrange(current_total), desired_number_tags);
	sample_list.sort();
	#print sample_list;
	count = 0;
	index = 0; 
	infile = open(raw_bed_file,'r');
	outfile = open(out_file_name, 'w');
	for line in infile:
		if count == sample_list[index]:	
			outfile.write(line);
			if index == desired_number_tags-1:
				break;
			else:
				index += 1;
		count +=1;
	outfile.close();
	infile.close()		


def main(argv):
	parser = OptionParser();
	parser.add_option("-f", "--rawtagfile", action="store", type="string",
			  dest="raw_bed_file", help="raw bed file",
			  metavar="<file>");
	parser.add_option("-n", "--desirednumberoftags", action="store", type="int",
			  dest="desired_number_tags", help="desired number of tags",
			  metavar="<int>");
    	parser.add_option("-o", "--slicedrawtagfile", action="store", type="string",
			  dest="out_file_name", help="sliced raw bed file",
			  metavar="<file>");
	(opt, args) = parser.parse_args(argv);
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	random_sample(opt.desired_number_tags, opt.raw_bed_file, opt.out_file_name);
	total = get_total_tag_counts.get_total_tag_counts(opt.out_file_name);
	print "The number of tags in " + opt.out_file_name + ' is ' + str(total);
	

if __name__ == "__main__":
	main(sys.argv)