#!/usr/bin/env python
# 
# Authors: Weiqun Peng
#
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

def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--islandsummaryfile", action="store", type="string",
			  dest="islandsummaryfile", metavar="<file>",
			  help="island summary file obtained from SICER")
	parser.add_option("-o", "--outfile", action="store", type="string",
			  dest="out_file", metavar="<file>",
			  help="island bed file")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	file = open(opt.islandsummaryfile,'r')
	ofile = open(opt.out_file, 'w')
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			ssline = sline[0:4]
			line =  '\t'.join(ssline) + '\n';
			ofile.write(line);
	file.close()
	ofile.close()

			
if __name__ == "__main__":
	main(sys.argv)
