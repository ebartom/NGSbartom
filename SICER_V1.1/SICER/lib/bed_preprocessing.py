#!/usr/bin/env python
# Authors: Weiqun Peng, Chongzhi Zang, Dustin E Schones and Keji Zhao
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

import BED;
import UCSC;
import bisect;
import GenomeData;
import Utility

Dir = os.getcwd();
grep = "/bin/grep";
cat = "/bin/cat";
plus = re.compile("\+");
minus = re.compile("\-");


def breakUpStrands(bed_list):
	"""
	input: a list of bed6 object 
	Return: two lists of bed6 objects, one for plus strand, one for minus strand. 
	"""

	plus_bed_list = [];
	minus_bed_list = [];
	for b in bed_list:
		if plus.match(b.strand):
			plus_bed_list.append(b);
		elif minus.match(b.strand):
			minus_bed_list.append(b);
	return (plus_bed_list, minus_bed_list)



def find_read_copy_distribution(sorted_bed_list):
	"""
	Input:  
		sorted_bed_list: a list of sorted bed6 objects. Already assumed that 
		the tags are from one chromosome and in one direction. 
	Return: the histogram of the tag copies
	"""
	assert (Utility.is_bed_sorted(sorted_bed_list) == 1)

	unique_tag_histogram = [0] * 100;
	if (len(sorted_bed_list) != 0):
		total_number_tags = len(sorted_bed_list);
		current_value = (sorted_bed_list[0]).start;
		current_count = 1;
		for index in range(1, len(sorted_bed_list)):
			item = sorted_bed_list[index];
			if (item.start != current_value):
				if (len(unique_tag_histogram)-1)<current_count:
					unique_tag_histogram +=[0]*(current_count-(len(unique_tag_histogram)-1));
				unique_tag_histogram[current_count] += 1;
				current_value = item.start;
				current_count = 1; #reset
			else:
				current_count += 1;
		#last read	
		if (len(unique_tag_histogram)-1)<current_count:
			unique_tag_histogram +=[0]*(current_count-(len(unique_tag_histogram)-1));
		unique_tag_histogram[current_count] += 1; 
	return unique_tag_histogram;

def find_n_copy_reads(sorted_bed_list, n):	
	"""
	Input:  
		sorted_bed_list: a list of sorted bed6 objects. Already assumed that 
				the tags are from one chromosome and in one direction. 
		n: the copies for a read 
	Return: the list of BED6 reads with copy number equal to n.
	"""
	
	assert (Utility.is_bed_sorted(sorted_bed_list) == 1)
	
	n_copy_read_list=[];
	temp_list = [];
	
	if (len(sorted_bed_list) != 0):
		total_number_tags = len(sorted_bed_list);
		temp_list.append(sorted_bed_list[0]);
		current_value = (sorted_bed_list[0]).start;
		current_count = 1;
		
		for index in range(1, len(sorted_bed_list)):
			item = sorted_bed_list[index];
			
			if (item.start != current_value):
				if (current_count==n):
					n_copy_read_list.extend(temp_list);
				current_value = item.start;
				current_count = 1; #reset
				temp_list = [];
				temp_list.append(item);
				
			else:
				current_count += 1;
				temp_list.append(item);
		#last read	
		if (current_count==threshold):
			n_copy_read_list.extend(temp_list);
	return n_copy_read_list;
			
def find_multi_copy_reads(sorted_bed_list, threshold):
	"""
	Input: 	
		sorted_bed_list: a list of sorted bed6 objects. Already assumed that 
						the tags are from one chromosome and in one direction. 
		threshold:	the threshold for read copy
	Return: the list of BED6 reads with copy number above or equal to threshold.
	"""
	multiple_copy_read_list=[];
	temp_list = [];
	
	assert (Utility.is_bed_sorted(sorted_bed_list) == 1)
	
	if (len(sorted_bed_list) != 0):
		total_number_tags = len(sorted_bed_list);
		current_value = (sorted_bed_list[0]).start;
		temp_list.append(sorted_bed_list[0]);
		current_count = 1;
		
		for index in range(1, len(sorted_bed_list)):
			item = sorted_bed_list[index];
			if (item.start != current_value):
				if (current_count>=threshold):
					#current_tag.score = current_count;
					#multiple_copy_read_list.append(current_tag);
					multiple_copy_read_list.extend(temp_list);
				current_value = item.start;
				current_count = 1; #reset
				temp_list = [];
				temp_list.append(item);
			else:
				current_count += 1;
				temp_list.append(item);
		#last read	
		if (current_count>=threshold):
			#item.score = current_count;
			#multiple_copy_read_list.append(item);
			multiple_copy_read_list.extend(temp_list);
	return multiple_copy_read_list;
	

def filter_reads(sorted_bed_list, cutoff, outfile):
	"""
	A read has n copies in the sorted_bed_list. If n<=cutoff, all the copies are retained.
	If n>cutoff, only cutoff number of copies of the read are retained.  
	
	Output: write bed objects with the extra redundant copies filtered out.If the number of reads in zero, then that file is not generated. 
	Return: the number of reads remained
	"""
	assert (Utility.is_bed_sorted(sorted_bed_list) == 1)
	counter2 = 0;
	if (len(sorted_bed_list) != 0):
		out = open(outfile, 'w')
		total_number_tags = len(sorted_bed_list);
		current_value = (sorted_bed_list[0]).start;
		current_count = 1;
		current_tag = sorted_bed_list[0];
		for index in range(1, len(sorted_bed_list)):
			item = sorted_bed_list[index];
			if (item.start != current_value):
				if (current_count <= cutoff):
					write(current_tag, out);
					counter2 +=1;
				current_value = item.start;
				current_count = 1; 
				current_tag = item;
			else:
				if (current_count <= cutoff):
					write(current_tag, out);
					counter2 +=1;
				current_count += 1;
				
		if (current_count <= cutoff): #last tag
			write(current_tag, out);
			counter2 +=1;	
		out.close();
	return counter2;


def filter_reads_add(sorted_bed_list, cutoff, outfile):
	"""
	A read has n copies in the sorted_bed_list. If n<=cutoff, all the copies are retained.
	If n>cutoff, only cutoff number of copies of the read are retained.  
	
	Output: write bed objects with the extra redundant copies filtered out. This function 
			uses out = open(outfile, 'a') instead of out = open(outfile, 'w')
	Return: the number of reads remained
	"""
	counter2 = 0;
	if (len(sorted_bed_list) != 0):
		out = open(outfile, 'a')
		total_number_tags = len(sorted_bed_list);
		current_value = (sorted_bed_list[0]).start;
		current_count = 1;
		current_tag = sorted_bed_list[0];
		for index in range(1, len(sorted_bed_list)):
			item = sorted_bed_list[index];
			if (item.start != current_value):
				if (current_count <= cutoff):
					write(current_tag, out);
					counter2 +=1;
				current_value = item.start;
				current_count = 1; 
				current_tag = item;
			else:
				if (current_count <= cutoff):
					write(current_tag, out);
					counter2 +=1;
				current_count += 1;
				
		if (current_count <= cutoff): #last tag
			write(current_tag, out);
			counter2 +=1;
		out.close();
	return counter2;

def write (item, out):
	"""
	write one line into outfile. The file openning and closing is handled by outside. 
	"""
	#chrom, start, end, name, score, strand
	outline = item.chrom + "\t" + str(item.start) + "\t" + str(item.end) + "\t" + item.name + "\t" + str(int(item.score)) + "\t" + item.strand + "\n";
	out.write(outline);
	
def write_list (bed6_list, out):
	"""
	write a bed6_list into outfile. The file openning and closing is handled from outside. 
	"""
	for item in bed6_list:
		#chrom, start, end, name, score, strand
		write (item, out);
	
	
def combine_histogram(a, b):
	t=[];
	if len(a)<len(b):
		t = b;
		for index in xrange(len(a)):
			t[index]  += a[index];
	else:
		t = a;
		for index in xrange(len(b)):
			t[index]  += b[index];
	return t;

def find_total_in_histogram(histogram, threshold=0):
	"""
	Threshold serves at the starting value for integration, inclusive.
	"""
	total = 0;
	for index in range(threshold, len(histogram)):
		total += index * histogram[index];
	return total;
	
def write_histogram(a, outfile):
	out = open(outfile, 'w');
	for index in xrange(len(a)):
		if (a[index] != 0):
			outline = str(index) + "\t" + str(a[index]) +"\n";
			out.write(outline); 
	out.close();
	
def combine_read_copy_distribution(species, file_name):
	"""
	file_name is for the raw tag file. 
	need BED6 to split the positive and negative tags. 
	"""
	chroms = GenomeData.species_chroms[species];
	histogram =[];
	bed_vals = BED.BED(species, file_name, "BED6", 0);
	for chrom in chroms:
		if chrom in bed_vals.keys():
			sorted_bed_list = (bed_vals[chrom]).sort(key=operator.attrgetter('start'));
			(plus_bed_list, minus_bed_list) = breakUpStrands(sorted_bed_list);
			plus_histogram = find_read_copy_distribution(plus_bed_list);
			histogram = combine_histogram(plus_histogram, histogram)
			minus_histogram = find_read_copy_distribution(minus_bed_list);
			histogram = combine_histogram(minus_histogram, histogram);
	return histogram;	

			
			

		
