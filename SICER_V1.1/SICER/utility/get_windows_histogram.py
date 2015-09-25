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

import GenomeData
import get_total_tag_counts

Dir = os.getcwd();

def get_total_num_windows (bed_graph_file, threshold = 0):
	"""
	Get total tag counts given the current experimental run
	file should be a summary bed graph file.
	can be used for bed summary or islands file
	
	return (the read counts on windows, and the number of windows) whose read counts >= threshold.
	"""
	count=0;
	infile = open(bed_graph_file,'r');
	for line in infile:
        	""" check to make sure not a header line """
        	if not re.match("track", line):
			line = line.strip();
            		sline = line.split();
			assert ( len(sline) == 4 );
			value = atof(sline[3]);
			if (value >= threshold):	
				count += 1;
	infile.close();
	return count;


def get_windows_histogram_from_bedsummary(species, summary_graph, window_size):
	"""
		Build the histogram for the windows according to the tag-count in each window. 
		It uses the summary.graph type that includes all the chromosomes. 
		Each line of the file should look like:
			chrom start end tag_count
	
	"""
	totalbp = 0.0;
	value_list = GenomeData.species_chrom_lengths[species].values();
	for values in value_list: totalbp += values;
	totalnumber_windows = round(totalbp/window_size);
	totalnumber_occupied_windows = 0.0;
	#print "Total number of windows is " + str(totalnumber_windows);
		
	current_max = 0;
	windows_hist = [0];
	
	for chrom in GenomeData.species_chroms[species]:
		for read in summary_graph[chrom]:
			read_count = int(read.value);
			if read_count > current_max:
				windows_hist += [0]*(read_count-current_max);
				current_max = read_count;
			windows_hist[read_count] += 1.0;
			totalnumber_occupied_windows += 1.0;	
	# The summary graph might or might not register those windows with zero tags.
	windows_hist[0] += totalnumber_windows - totalnumber_occupied_windows;
	
	#Get total number of tags
	total_number_tags = 0;
	for i in range (0, len(windows_hist)):
		total_number_tags += i*windows_hist[i];
	print "Total number of tags is " + str(total_number_tags);
	return windows_hist;

def output_windows_histogram(windows_hist, output_filename):
	total_num_windows = sum(windows_hist);
	print "Total number of windows is " + str(total_num_windows);	
	if (output_filename != ""):
		outfile = open(output_filename, 'w');
		outline = "# read_cout"+ "\t" +"Observed_histogram" + "\n";
		outfile.write(outline);
		for i in range (0, len(windows_hist)):
			outline = str(i) + "\t" + str(windows_hist[i]) +  "\n";
			outfile.write(outline);
		outfile.close();

def get_windows_histogram(species, summary_graph_file, window_size, rescale_factor=1):

	"""
	The rescale factor is used to rescale the total tag count.
	Build the histogram for the windows according to the tag-count in each window. 
	It uses the summary.graph type file  that includes all the chromosomes. 
	Each line of the file should look like:
	chrom start end tag_count
	"""
	totalbp = 0.0;
	value_list = GenomeData.species_chrom_lengths[species].values();
	for values in value_list: totalbp += values;
	totalnumber_windows = round(totalbp/window_size);
	totalnumber_occupied_windows = 0.0;
	#print "Total number of windows is " + str(totalnumber_windows);
		
	current_max = 0;
	windows_hist = [0];

	infile = open(summary_graph_file, 'r');
		
        for line in infile:
        	""" check to make sure not a header line """
		if not re.match("track", line):
			line = line.strip();
			sline = line.split();
			assert (len(sline) == 4);
			tags_in_window = atof(sline[3]);
			tags_in_window = int(rescale_factor*tags_in_window)
			if tags_in_window > current_max:
				windows_hist += [0]*(tags_in_window-current_max);
				current_max = tags_in_window;
			windows_hist[tags_in_window] += 1.0;
			totalnumber_occupied_windows += 1.0;
	# The summary graph might or might not register those windows with zero tags.		
	windows_hist[0] += totalnumber_windows - totalnumber_occupied_windows;

	#Get total number of tags
	total_number_tags = 0;
	for i in range (0, len(windows_hist)):
		total_number_tags += i*windows_hist[i];
	print "Total number of tags is " + str(total_number_tags);

	return windows_hist;
	

    
def main(argv):
	
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="species under consideration", metavar="<str>")
	parser.add_option("-b", "--_summary_graph_file", action="store", type="string",
                      dest="summary_graph_file", help="summary graph file of", metavar="<file>")
	parser.add_option("-w", "--window_size", action="store", type="int",
                      dest="window_size", help="window size", metavar="<int>")
	parser.add_option("-o", "--output_filename", action="store", type="string",
                      dest="output_filename", help="output_filename", metavar="<file>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
  	if opt.species not in GenomeData.species_chroms.keys():
		print "The species is not recognized!!";
	else:
		hist = get_windows_histogram(opt.species, opt.summary_graph_file, opt.window_size);
		output_windows_histogram(hist, opt.output_filename);
    

if __name__ == "__main__":
	main(sys.argv)


        
