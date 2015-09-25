#!/usr/bin/env python
# Copyright (c) 2010 The George Washington University
# Copyright (c) 2007 NHLBI, NIH
# Authors: Weiqun Peng, Chongzhi Zang and Keji Zhao
#
# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
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

def get_tag_counts_by_chrom(tag_bed_file, chrom):
	"""
	Get total tag counts given the current experimental run
	file should be a bed file.
	
	chrom: string
	"""
	counts =0.0;
	infile = open(tag_bed_file,'r');
	for line in infile:
        	""" check to make sure not a header line """
        	if not re.match("track", line):
			line = line.strip();
            		sline = line.split();
			if (sline[0] == chrom):
				counts += 1.0;
	infile.close();
	return counts;



def get_total_tag_counts(tag_bed_file):
	"""
	Get total tag counts given the current experimental run
	file should be a bed file.
	"""
	counts =0.0;
	infile = open(tag_bed_file,'r');
	for line in infile:
        	""" check to make sure not a header line """
        	if not re.match("track", line):
				counts += 1.0;
	infile.close();
	return counts;

def get_total_tag_counts_bed_graph(summary_graph_file, bed_val={}, threshold = 0):
	"""
	Get total tag counts given the current experimental run
	file should be a summary bed graph file.
	can be used for bed summary or islands file
	
	"""
	count = 0.0;
	if summary_graph_file != "":
		infile = open(summary_graph_file,'r');
		for line in infile:
			""" check to make sure not a header line """
			if not re.match("track", line):
				line = line.strip();
				sline = line.split();
				assert ( len(sline) == 4 );
				value = atof(sline[3]);
				if (value >= threshold):	
					count += value;
		infile.close();
	elif len(bed_val)>0:
		for chrom in bed_val.keys():
			for item in bed_val[chrom]:
				if (item.value >= threshold):
					count +=item.value;
	else:
		print "wrong input!!"; 
	return count;


def main(argv):
	parser = OptionParser();
	parser.add_option("-f", "--tagfile", action="store", type="string",
			  dest="tagfile", help="file with tag coords in bed format",
			  metavar="<file>");
    
	(opt, args) = parser.parse_args(argv);
	if len(argv) < 2:
		parser.print_help()
		sys.exit(1)

	total = get_total_tag_counts(opt.tagfile);    	
	print "The number of tags in " + opt.tagfile + ' is ' + str(total);
	

if __name__ == "__main__":
	main(sys.argv)
