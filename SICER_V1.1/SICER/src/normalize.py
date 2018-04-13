#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


def normalize_tag_count(input_file, ColumnIndex, total, output_file):
	"""
	Column index is 0 based 
	"""
	infile = open(input_file,'r')
	outfile = open(output_file, 'w')
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) >= (ColumnIndex + 1):
				sline[ColumnIndex] = str(atof(sline[ColumnIndex])/total)
				outline = '\t'.join(sline) +'\n'
				outfile.write(outline)
	infile.close()
	outfile.close()


def total_counts(file, column):
	"""
	Column index is 0 based 
	"""
	total = 0.0
	infile = open(file,'r')
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) >= (column + 1):
				total += atof(sline[column])
	infile.close()
	return total


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original data file to be normalized")
	parser.add_option("-a", "--tagcountcolumn", action="store", type="int", dest="colum", metavar="<int>", help="colum index for tag counts to be normalized, start from 0")
	parser.add_option("-t", "--ScalingFactor", action="store", type="float", dest="ScalingFactor", metavar="<float>", help="Scaling factor for normalization, for example, scaling factor is 1000000 if normalization by total read count per million")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	total = total_counts(opt.input_file, opt.colum) / opt.ScalingFactor
	normalize_tag_count(opt.input_file, opt.colum, total, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)
