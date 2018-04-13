#!/usr/bin/env python
# 
# Authors: Chongzhi Zang, Weiqun Peng
#
# This module compares two libraries on islands, which by definition are NON-OVERLAPPING! It can not be used for comparison of 
# libraries on genes, which are potentially overlapping!

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

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import BED
import GenomeData;
import associate_tags_with_regions
import SeparateByChrom
import get_total_tag_counts
import Utility
import scipy.stats

def pvaule (chip_read_count, control_read_count, scaling_factor, pseudo_count):
	"""
	Currently using poisson distribution
	
	scaling_factor: the factor that accounts for the differences of control library and ChIP library. effective control read count
	is control_read_count * scaling factor
	pseudocount: when control_read_count is zero, replace zero with pseudocount to alleviate the impact of statistical fluctuation
	
	output: pvalue
	"""
	if control_read_count > 0:
		average = control_read_count * scaling_factor;
	else:
		average = pseudo_count * scaling_factor
	if  chip_read_count > average:
		pvalue = scipy.stats.poisson.sf(chip_read_count, average)[()]; 
	else:
		pvalue = 1;
	return pvalue;
	
def fdr(pvalue_list):
	"""
	Calculate the multiple testing corrected p-value using BH
	"""	
	fdr_list=[];
	pvaluearray=scipy.array(pvalue_list);
	totalnumber = len(pvalue_list);
	pvaluerankarray=scipy.stats.rankdata(pvaluearray);
	for i in range(totalnumber):
		fdr_value = pvalue_list[i] * totalnumber/pvaluerankarray[i];
		if fdr_value > 1 :
			fdr_value = 1
		fdr_list.append(fdr_value);
	return fdr_list;
		

def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18, etc", metavar="<str>")
	parser.add_option("-a", "--rawreadfileA", action="store", type="string", dest="readfileA", metavar="<file>", help="raw read file A in bed format")
	parser.add_option("-b", "--rawreadfileB", action="store", type="string", dest="readfileB", metavar="<file>", help="raw read file B in bed format")
	parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", metavar="<int>", help="average size of a fragment after A experiment")
	parser.add_option("-d", "--islandfile", action="store", type="string", dest="islandfile", metavar="<file>", help="island file in BED format")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="island read count summary file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
        	parser.print_help()
        	sys.exit(1)
		
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	if not Utility.fileExists(opt.readfileA):
		print opt.readfileA, " not found";
		sys.exit(1)
	if not Utility.fileExists(opt.readfileB):
		print opt.readfileB, " not found";
		sys.exit(1)	
	
	A_library_size = get_total_tag_counts.get_total_tag_counts(opt.readfileA);
	B_library_size = get_total_tag_counts.get_total_tag_counts(opt.readfileB);
	print "Library size of ", opt.readfileA, ":  ", A_library_size
	print "Library size of ", opt.readfileB, ":  ", B_library_size
	
	totalA = 0;
	totalB = 0;
	
	islands = BED.BED(opt.species, opt.islandfile, "BED3", 0);
	
	# separate by chrom the A library
	SeparateByChrom.separateByChrom(chroms, opt.readfileA, '.bed1');
	# separate by chrom the B library
	SeparateByChrom.separateByChrom(chroms, opt.readfileB, '.bed2');
	
	
	island_A_readcount = {};
	island_B_readcount = {};
	
	#Find read counts on the islands
	for chrom in chroms:
		if chrom in islands.keys():
			if len(islands[chrom]) != 0:
				island_list = islands[chrom];
				if Utility.is_bed_sorted(island_list) == 0:
					island_list.sort(key=operator.attrgetter('start'));
					
				island_start_list = []
				island_end_list = []
				for item in island_list:
					island_start_list.append(item.start)
					island_end_list.append(item.end)
	
				island_A_readcount_list=[0]*len(island_list);
				read_file = chrom + ".bed1";
				f = open(read_file,'r')
				for line in f:
					if not re.match("#", line):
						line = line.strip()
						sline = line.split()
						position = associate_tags_with_regions.tag_position(sline, opt.fragment_size)
						index =associate_tags_with_regions.find_readcount_on_islands(island_start_list, island_end_list, position);
						if index >= 0:
							island_A_readcount_list[index] += 1;
							totalA += 1;
				f.close();
				island_A_readcount[chrom] = island_A_readcount_list;
							
				island_B_readcount_list=[0]*len(island_list);
				read_file = chrom + ".bed2";
				f = open(read_file,'r')
				for line in f:
					if not re.match("#", line):
						line = line.strip()
						sline = line.split()
						position = associate_tags_with_regions.tag_position(sline, opt.fragment_size)
						index = associate_tags_with_regions.find_readcount_on_islands(island_start_list, island_end_list, position);
						if index >= 0:
							island_B_readcount_list[index] += 1;
							totalB += 1;
				f.close();		
				island_B_readcount[chrom] = island_B_readcount_list;			
						
	#A_background_read = A_library_size - totalA;
	#B_background_read = B_library_size - totalB;
	
	print "Total number of A reads on islands is: ", totalA; 
	print "Total number of B reads on islands is: ", totalB; 

	# Calculate the p value.
	library_scaling_factor = A_library_size*1.0/B_library_size; #A vs B
	pseudo_count = 1; 
	pvalue_A_vs_B_list = [];
	pvalue_B_vs_A_list = [];
	for chrom in chroms:
		if chrom in islands.keys():
			if len(islands[chrom]) != 0:
				island_list = islands[chrom];
				for index in xrange(len(island_list)):
					item = island_list[index];
					Acount = (island_A_readcount[chrom])[index]; 
					Bcount = (island_B_readcount[chrom])[index];
					pvalue_A_vs_B = pvaule (Acount, Bcount, library_scaling_factor, pseudo_count);
					pvalue_A_vs_B_list.append(pvalue_A_vs_B);
					pvalue_B_vs_A = pvaule (Bcount, Acount, 1/library_scaling_factor, pseudo_count);
					pvalue_B_vs_A_list.append(pvalue_B_vs_A);
	#Calculate the FDR
	fdr_A_vs_B_list = fdr(pvalue_A_vs_B_list);
	fdr_B_vs_A_list = fdr(pvalue_B_vs_A_list);


	#Output the islands read counts, normalized read counts, fc, pvalue both ways
	scaling_factor = 1000000; 
	out = open(opt.out_file, 'w');
	outline = '#chrom' + "\t" + 'start' + "\t" + 'end' + "\t" + "Readcount_A" + "\t" + 'Normalized_Readcount_A' + "\t" + 'ReadcountB' + "\t" + 'Normalized_Readcount_B' + "\t" + "Fc_A_vs_B" + "\t" + "pvalue_A_vs_B" + "\t" + "FDR_A_vs_B" + "\t" + "Fc_B_vs_A" + "\t" + "pvalue_B_vs_A" + "\t" + "FDR_B_vs_A"  + "\n"; 	
	out.write(outline);
	ii=0;
	for chrom in chroms:
		if chrom in islands.keys():
			if len(islands[chrom]) != 0:
				island_list = islands[chrom];
				for index in xrange(len(island_list)):
					item = island_list[index];
					Acount = (island_A_readcount[chrom])[index]; 
					Bcount = (island_B_readcount[chrom])[index];
					normalized_A = Acount/ float(A_library_size) * scaling_factor;
					normalized_B = Bcount/ float(B_library_size) * scaling_factor;
					fc_A_vs_B = ((Acount + pseudo_count)*1.0/(Bcount + pseudo_count))/library_scaling_factor;
					fc_B_vs_A = ((Bcount + pseudo_count)*1.0/(Acount + pseudo_count)) * library_scaling_factor;
					outline = item.chrom + "\t" + str(item.start) + "\t" + str(item.end) + "\t" + str(Acount) + "\t"  +  str(normalized_A) + "\t"  +  str(Bcount) + "\t" + str(normalized_B) + "\t" +  str(fc_A_vs_B) + "\t" + str(pvalue_A_vs_B_list[ii]) + "\t" + str(fdr_A_vs_B_list[ii]) + "\t" + str(fc_B_vs_A) + "\t" + str(pvalue_B_vs_A_list[ii]) + "\t" + str(fdr_B_vs_A_list[ii]) + "\n";	
					out.write(outline);
					ii += 1;		
	out.close();

	SeparateByChrom.cleanup(chroms, '.bed1');
	SeparateByChrom.cleanup(chroms, '.bed2');


	# Calculate the correlations using normalized read counts
	A_array=();
	B_array=();
	for chrom in chroms:
		if chrom in islands.keys():
			if len(islands[chrom]) != 0:
				temp_array= scipy.array(island_A_readcount[chrom]);
				A_array=scipy.concatenate((temp_array, A_array));
				temp_array= scipy.array(island_B_readcount[chrom]);
				B_array=scipy.concatenate((temp_array, B_array));
	#Normalization to reads per million
	A_array = A_array/float(A_library_size) * scaling_factor;
	B_array = B_array/float(B_library_size) * scaling_factor;
	pearson=scipy.stats.pearsonr(A_array, B_array);
	print "Pearson's correlation is: ", pearson[0], " with p-value ",  pearson[1];
	spearman = scipy.stats.spearmanr(A_array, B_array);
	print "Spearman's correlation is: ", spearman[0], " with p-value ",  spearman[1];
	
	



if __name__ == "__main__":
	main(sys.argv)
