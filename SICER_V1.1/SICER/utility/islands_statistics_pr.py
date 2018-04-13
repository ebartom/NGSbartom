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

import get_total_tag_counts
import GenomeData
import filter_summary_graphs
import BED

def output_histogram(bins, histogram, filename):
	outfile = open(filename, 'w');
	assert (len(bins) == len(histogram));
	for i in range(0, len(histogram)):
		if histogram[i]>0:
			outline = str(bins[i]) + "\t " + str(histogram[i]) + "\n";
			outfile.write(outline);
	outfile.close();
	
	
def find_islands_length_histogram(bed_vals, islands_file, outfilename):
	"""
	For use of obtaining the length histogram of islands. 
	"""
	current_max = 0;
	total_number_islands = 0.0;
	total_length = 0;
	histogram = [0];

	if ( bed_vals != {} and islands_file == ""):
		for chrom in bed_vals.keys():
			for bed_item in bed_vals[chrom]:
				length = bed_item.end - bed_item.start + 1;
				assert (length >= 0);
				total_length += length;
				if length > current_max:
					histogram += [0]*(length-current_max + 1);
					current_max = len(histogram);
				histogram[length] += 1.0;
				total_number_islands += 1.0;
	elif ( bed_vals == {} and islands_file != ""):
		infile = open(islands_file, 'r');
		for line in infile:
			if not re.match("track", line):
				line = line.strip();
				sline = line.split();
				assert (len(sline) == 4);
				length = atoi(sline[2]) - atoi(sline[1]) + 1;
				assert (length >= 0);
				total_length += length;
				if length > current_max:
					histogram += [0]*(length-current_max + 1);
					current_max = len(histogram);
				histogram[length] += 1.0;
				total_number_islands += 1.0;
		infile.close();
	else:
		print "wrong input!";	
	
	outfile = open(outfilename, 'w');
	outline = "# The totoal length of the islands is: " + str(total_length) + "\n";
	outfile.write(outline);
	
	normalization =0.0;
	for i in range(0, len(histogram)):
		if histogram[i]>0:
			normalization += (float(histogram[i])/total_number_islands);
			outline = str(i) + "\t " + str(histogram[i]) + "\n";
			outfile.write(outline);
	outfile.close();
	if (total_number_islands >0): assert (fabs(normalization-1)<.00000000001);
	return total_length;

	
def find_islands_score_histogram(bed_vals, islands_file, bin_size, outfile):
	"""
	No normalization is applied
	The input can be either bed_vals a BED object or a file.
	"""	
	total_number_islands = 0.0;
	total_score = 0.0; 
	histogram = [0];
	
	if ( bed_vals != {} and islands_file == ""):
		for chrom in bed_vals.keys():
			for item in bed_vals[chrom]:
				index = int(item.value/bin_size);
				if index >= len(histogram):
					histogram += [0] * (index-len(histogram)+1)
				histogram[index] +=1.0;
				total_score += item.value;
				total_number_islands += 1.0;
	elif ( bed_vals == {} and islands_file != ""):
		infile = open(islands_file, 'r');
		for line in infile:
			""" check to make sure not a header line """
			if not re.match("track", line):
				line = line.strip();
				sline = line.split();
				assert (len(sline) == 4);
				score = atof(sline[3]);
				index = int(score/bin_size);
				if index >= len(histogram):
					histogram += [0]*(index-len(histogram)+1);
				histogram[index] += 1.0;
				total_score += score;
				total_number_islands += 1.0;
		infile.close();
	else:
		print "wrong input!!!!!!!!!!!!!!!!!!!!!";	
	# print "Total number of islands is :",  total_number_islands;
	
	bins=[];
	for index in xrange(len(histogram)):
		# [0, bin_size) in first bin, [bin_size, 2*bin_size) in second bin, etc
		bins.append(bin_size*index);
	if outfile != "":
		output_histogram(bins, histogram, outfile);
		
	return total_score;

def  get_island_read_counts (species, summary_graph_file, islands_file,  window_size, out_file, window_read_count_threshold=0):
	"""
	Filter summary graphs using the islands. Find the read count on islands
	"""
	total_read_count = 0;
	
	assert (species in GenomeData.species_chroms.keys())
	windows_on_islands = filter_summary_graphs.find_windows_on_islands(species, summary_graph_file, islands_file,  window_size, out_file, window_read_count_threshold);
	
	total_read_count = get_total_tag_counts.get_total_tag_counts_bed_graph("", windows_on_islands, window_read_count_threshold);
	return total_read_count;
	
	
def windowlize_islands(species, islands_file, window_size, shift):
	"""
	For use in dead zone analysis.
	"""
	windows = {}; 
	infile = open(islands_file, 'r');
	for line in infile:
		if not re.match("track", line):
			line = line.strip();
			sline = line.split();
			chrom_length = ((GenomeData.species_chrom_lengths)[species])[sline[0]];
			
			if sline[0] not in windows.keys():
				num_windows = int(ceil(float(chrom_length)/window_size));
				windows[sline[0]]=[0]*num_windows;
			
			# start and end are both from 5".
			start = atoi(sline[1]) + shift;
			end = atoi(sline[2]) - 1 + shift; # BED entires are half-open intervals.
			if (end >= 0) and (start <= (chrom_length-1)): 
				if start < 0:
					start = 0;
				if end > chrom_length-1:
					end =  chrom_length-1;
				start_window = start/window_size;
				end_window = end/window_size;
				assert (end_window>=start_window);
				if start_window == end_window:
					windows[sline[0]][start_window] += end - start + 1;
				elif start_window < end_window:
					for index in range(start_window, end_window+1):		
						if index == start_window:
							windows[sline[0]][index] += window_size * (index+1) - start;
						elif index == end_window:
							windows[sline[0]][index] += end - index*window_size + 1;
						else:
							windows[sline[0]][index] += window_size;
	return windows;	 						
						

	
	
def main(argv):
	
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="mm8, hg18, background, etc", metavar="<str>")
	parser.add_option("-b", "--summary_graph_file", action="store", type="string",
                      dest="summary_graph_file", help="bed summary graph file of", metavar="<file>")
	parser.add_option("-i", "--islands_file", action="store", type="string",
                      dest="islands_file", help="islands file", metavar="<file>")
	parser.add_option("-w", "--window_size", action="store", type="int",
                      dest="window_size", help="window size of summary", metavar="<int>")     
	parser.add_option("-o", "--islands_score_histogram_file", action="store", type="string",
                      dest="islands_score_histogram_file", help="islands histogram file", metavar="<file>")
	parser.add_option("-q", "--islands_length_histogram_file", action="store", type="string",
                      dest="islands_length_histogram_file", help="islands length histogram file", metavar="<file>")
	parser.add_option("-r", "--island_filtered_summary_graph", action="store", type="string",
                      dest="island_filtered_summary_graph", default = "", help=" Optional. The default is not to do it. ", metavar="<file>")	   	   

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
        	parser.print_help()
        	sys.exit(1)
	if opt.species in GenomeData.species_chroms.keys():
		genome_length = sum ( GenomeData.species_chrom_lengths[opt.species].values());
		chroms = GenomeData.species_chroms[opt.species];
		
		total_tag_counts = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.summary_graph_file);
		print "Total read count is:" , total_tag_counts;
		
		total_length = find_islands_length_histogram({}, opt.islands_file, opt.islands_length_histogram_file);
		print "Total islands length is: ", total_length, ";      Length coverage = total_length_of_islands/genome_length is: ", total_length*1.0/genome_length;
		
		bin_size=0.1;
		total_score = find_islands_score_histogram({}, opt.islands_file, bin_size, opt.islands_score_histogram_file);
		print "Total islands score is: ", total_score;
		
		if (opt.island_filtered_summary_graph != ""):
			read_count_on_islands = get_island_read_counts (opt.species, opt.summary_graph_file, opt.islands_file,  opt.window_size, opt.island_filtered_summary_graph, 0)
			print "Total read count on island is: ", read_count_on_islands, " Read count coverage=read_count_on_islands/Total-read-count: ", read_count_on_islands/float(total_tag_counts); 
		
    	else: 
		print "This species is not in my list!"; 

if __name__ == "__main__":
	main(sys.argv)

