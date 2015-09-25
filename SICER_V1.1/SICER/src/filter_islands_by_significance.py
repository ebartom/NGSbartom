#!/usr/bin/env python
# 
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
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--islandsummary", action="store", type="string", dest="islandsummary", metavar="<file>", help="island summary file")
	parser.add_option("-p", "--significance", action="store", type="float", dest="significance", metavar="<float>", help="significance, default value -1", default=-1)
	parser.add_option("-c", "--columnindex", action="store", type="int", dest="column", metavar="<int>", help="column index for significance")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="File storing significant islands given significance")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
		
	inputfile = open(opt.islandsummary,'r');
	outfile = open(opt.out_file, 'w');

	totalislands = 0;	

	for line in inputfile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if atof(sline[opt.column]) <= opt.significance:
				totalislands += 1;
				outfile.write('\t'.join(sline)+'\n');
	print "Given significance",  opt.significance, ",  there are", totalislands, "significant islands";
	inputfile.close()
	outfile.close()
	
	
	
if __name__ == "__main__":
	main(sys.argv)