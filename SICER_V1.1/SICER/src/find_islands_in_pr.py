#!/usr/bin/env python
# Authors: Chongzhi Zang, Weiqun Peng, Dustin E Schones and Keji Zhao
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
# Version 1.1 11/9/2010


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import BED
import GenomeData
import get_total_tag_counts
import Background_island_probscore_statistics
import Utility

""" 
Take in coords for bed_gaph type summary files and find 'islands' of modifications.
There are a number options here that can be turned on or off depending on need
Right now:

(1) scan for all 'islands' where the single window or consecutive
windows

(2) remove all single window 'islands' -- this can be commented out
when looking for more localized signals (such as TF binding sites?)

(3) try to combine islands that are within gap distance of each other.
This gap distance is supplemented by a window_buffer just so we don't
do anything stupid with window sizes

(4) Remove all single window combined islands -- if step (2) was done,
this is redundant

(5) Lastly, filter out all the islands we've found that have a total
score < islands_minimum_tags
"""


# Factorial
def fact(m):
	value = 1.0;
	if m != 0:
		while m != 1:
			value = value*m;
			m = m - 1;
	return value;

# Return the log of a factorial, using Srinivasa Ramanujan's approximation when m>=20
def factln(m):
	if m<20:  
		value = 1.0;
		if m != 0:
			while m != 1:
				value = value*m;
				m = m - 1;
		return log(value);
	else:
		return m*log(m) -m + log(m*(1+4*m*(1+2*m)))/6.0 + log(pi)/2;

def poisson(i, average):
	if i<20:
		return exp(-average) * average**i / fact(i);
	else:
		exponent = -average + i*log(average) - factln(i);
		return exp(exponent);
	

def find_threshold(pvalue, average):
	"""
	Returns the thershold value T given the p-value requirement.
	Namely, P(T)+P(T+1)+ ... +P(infty) <p-vlue
	"""
	value = 1;
	index = 0;
	value -= poisson(index, average);
	while value > pvalue:
		index += 1;
		value -= poisson(index, average);	
	# index I is the value of which 1-P(0)-P(1)-...-P(I) < pvalue
	return index+1;



def removeSingleWindowIslands(islands, window):
    filtered_islands = [];
    for i in islands:
        size = i.end - i.start;
        if size > window:
            filtered_islands.append(i);
    return filtered_islands;



def combine_proximal_islands(islands, gap, window_size_buffer=3):
	"""
	islands: a list of BEd_GRAPH object: (chrom, start, end, value)
	
	Extend the regions found in the find_continuous_region function.
	If gap is not allowed, gap = 0, if one window is allowed, gap = window_size (200) 
	
	Return a list of combined regions.
	"""
	#print len(islands);
	proximal_island_dist = gap + window_size_buffer;
	Final_islands=[]
	if len(islands) > 0:
		if not Utility.is_bed_sorted (islands): islands.sort(key=operator.attrgetter('start'));
		current_island=islands[0];
		#print current_island;
		if len(islands) == 1:
			Final_islands = islands;
		else:
			for index in range(1, len(islands)):
				dist = islands[index].start - current_island.end;
				if dist <= proximal_island_dist:
					current_island.end = islands[index].end;
					current_island.value += islands[index].value;
				else:
					Final_islands.append (current_island);
					current_island = islands[index];
			# The last island:
			Final_islands.append (current_island);
	#print len(Final_islands);
	return Final_islands;


def find_region_above_threshold(island_list, islands_minimum_tags):
    filtered_islands = [];
    for island in island_list:
        if island.value >= (islands_minimum_tags-.0000000001): filtered_islands.append(island);
    return filtered_islands;


def find_region_above_threshold_from_file(Islands_file, species, islands_minimum_tags, out_islands_file):
    bed_vals = BED.BED(species, Islands_file, "BED_GRAPH");
    outputfile = open(out_islands_file, 'w');
    total_number_islands = 0.0;
    total_tags_on_islands = 0.0;
    for chrom in bed_vals.keys():
        islands = bed_vals[chrom];
        islands = find_region_above_threshold(islands, islands_minimum_tags);
        total_number_islands += len(islands);
        for i in islands:
            outline = chrom + " " + str(i.start) + " " + str(i.end) + " " \
                      + str(i.value) + "\n";	
            outputfile.write(outline);
            total_tags_on_islands += i.value;
    outputfile.close();



def main(argv):
	"""
	Probability scoring with random background model.
	
	"""
	parser = OptionParser()
	
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="mm8, hg18, background, etc", metavar="<str>")
	parser.add_option("-b", "--summarygraph", action="store",type="string", dest="summarygraph", help="summarygraph", metavar="<file>")
	parser.add_option("-w", "--window_size(bp)", action="store", type="int", dest="window_size", help="window_size(in bps)", metavar="<int>")
	parser.add_option("-g", "--gap_size(bp)", action="store", type="int",  dest="gap", help="gap size (in bps)", metavar="<int>")
	parser.add_option("-t", "--mappable_fraction_of_genome_size ", action="store", type="float", dest="fraction", help="mapable fraction of genome size", metavar="<float>")
	parser.add_option("-e", "--evalue ", action="store", type="float", dest="evalue", help="evalue that determines score threshold for significant islands", metavar="<float>")
	parser.add_option("-f", "--out_island_file", action="store", type="string", dest="out_island_file", help="output island file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 14:
        	parser.print_help()
        	sys.exit(1)

	if opt.species in GenomeData.species_chroms.keys():
		print "Species: ", opt.species;
		print "Window_size: ", opt.window_size;
		print "Gap size: ", opt.gap;
		print "E value is:", opt.evalue;
		
		total_read_count = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.summarygraph);
		print "Total read count:", total_read_count
		genome_length = sum (GenomeData.species_chrom_lengths[opt.species].values());
		print "Genome Length: ", genome_length;
		genome_length = int(opt.fraction * genome_length);

		average = float(total_read_count) * opt.window_size/genome_length; 
		print "Effective genome Length: ", genome_length;
		print "Window average:", average;
		
		window_pvalue = 0.20;
		bin_size = 0.001;
		print "Window pvalue:", window_pvalue;
		background = Background_island_probscore_statistics.Background_island_probscore_statistics(total_read_count, opt.window_size, opt.gap, window_pvalue, genome_length, bin_size);
		min_tags_in_window = background.min_tags_in_window
		print "Minimum num of tags in a qualified window: ", min_tags_in_window
		
		print "Generate the enriched probscore summary graph and filter the summary graph to get rid of ineligible windows "; 
		#read in the summary graph file
		bed_val = BED.BED(opt.species, opt.summarygraph, "BED_GRAPH");
		
		#generate the probscore summary graph file, only care about enrichment
		#filter the summary graph to get rid of windows whose scores are less than window_score_threshold
		
		filtered_bed_val = {};
		
		for chrom in bed_val.keys():
			if len(bed_val[chrom])>0:
				filtered_bed_val [chrom]= [];
				for index in xrange(len(bed_val[chrom])):
					read_count = bed_val[chrom][index].value;
					if ( read_count < min_tags_in_window):
						score = -1;
						#score = 0;
					else:
						prob = poisson(read_count, average);
						if prob <1e-250:
							score = 1000; #outside of the scale, take an arbitrary number.
						else:
							score = -log(prob);
					bed_val[chrom][index].value = score;
					if score > 0:
						filtered_bed_val[chrom].append( (bed_val[chrom])[index] );
					#print chrom, start, read_count, score;
		
		#write the probscore summary graph file
		#Background_simulation_pr.output_bedgraph(bed_val, opt.out_sgraph_file);
		
		#Background_simulation_pr.output_bedgraph(filtered_bed_val, opt.out_sgraph_file+".filtered");
		
		print "Determine the score threshold from random background"; 
		#determine threshold from random background
		hist_outfile="L" + str(genome_length) + "_W" +str(opt.window_size) + "_G" +str(opt.gap) +  "_s" +str(min_tags_in_window) + "_T"+ str(total_read_count) + "_B" + str(bin_size) +"_calculatedprobscoreisland.hist";
		score_threshold = background.find_island_threshold(opt.evalue); 
		# background.output_distribution(hist_outfile);
		print "The score threshold is: ", score_threshold;
		
		
		print "Make and write islands";
		total_number_islands = 0;
		outputfile = open(opt.out_island_file, 'w');
		for chrom in filtered_bed_val.keys():
			if len(filtered_bed_val[chrom])>0:
				islands = combine_proximal_islands(filtered_bed_val[chrom], opt.gap, 2);
				islands = find_region_above_threshold(islands, score_threshold);
				total_number_islands += len(islands);
				if len(islands)>0:
					for i in islands:
						outline = chrom + "\t" + str(i.start) + "\t" + str(i.end) + "\t" + str(i.value) + "\n";	
						outputfile.write(outline);
				else:
					print "\t", chrom, "does not have any islands meeting the required significance";
		outputfile.close();	
		print "Total number of islands: ", total_number_islands;
		
	else:
		print "This species is not in my list!"; 

if __name__ == "__main__":
	main(sys.argv)
