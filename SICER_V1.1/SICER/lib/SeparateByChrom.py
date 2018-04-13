#!/usr/bin/env python
# Copyright (c) 2007 NHLBI, NIH
# Authors: Dustin E Schones and Keji Zhao
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
import operator

import Utility

grep = "grep";
cat = "cat";
#grep = "/bin/grep";
#cat = "/bin/cat";

plus = re.compile("\+");
minus = re.compile("\-");


def separateByChrom(chroms, file, extension):
	"""
	It is ok if the chroms do not include all those existing in the file.
	"""
    	for chrom in chroms:
        	match = chrom + "[[:space:]]";
        	tmpFile = chrom + extension;
        	try:
            		if os.system('%s %s %s > %s' %
                         	(grep, match, file, tmpFile)): raise
        	except: sys.stderr.write( "Warning: " + str(chrom) + " reads do not exist in " + str(file) + "\n");


def combineAllGraphFiles(chroms, extension, final_out):
	"""
		Combine the seperately processed chromosomes, return the output file name
	"""
    	outfile = open(final_out,'w');
    	outfile.close();
    
    	for chrom in chroms:
		file = chrom + extension;
		if Utility.fileExists(file):
			try:
				if os.system('%s %s >> %s' %
					(cat, file, final_out)): raise
			except: 
				sys.stderr.write( "")
	#			sys.stderr.write("cat failed\n")
		else:
			print file, " file does not exist."
	return final_out


def cleanup(chroms, extension):
    for chrom in chroms:
        file = chrom + extension;
        try:
            if os.remove('%s' %
                         (file)): raise
        except: 
		    sys.stderr.write("")
#		    sys.stderr.write("clean up failed\n");

