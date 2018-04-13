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
# schonesde@mail.nih.gov).


"""
This module contains all the classes to deal with BED files
"""
import re, os, shutil, time, sys
from math import *
from string import *
import GenomeData;


plus = re.compile('\+');
minus = re.compile('\-');

bedError = "Error in BED class";

"""
  BED types:

  BED2: start, strand
  BED3: chrom, start and end
  BED6: BED3 + name + score + strand ('+' or '-')
  BED_GRAPH: bed graph format to mimic wiggle format:
               chrom, start, end, value
"""


#----------------------------
#----------------------------

class BED2:
    """
    Class for bed lines with 2 values: start and strand, this will be
    useful for things like TSS information or tags where the only
    important information is the start and strand.
    """
    def __init__(self, start, strand):
        self.start = start;
        self.strand = strand;

    def __set__(self, start, strand):
        self.start = start;
        self.strand = strand;

    def getCoord(self):
        outstring = str(self.start) + " " + self.strand;
        try:
            return outstring;
        except:
            sys.stderr.write("No coord information for %s\n" % self)
            return ''


#----------------------------
#----------------------------

class BED3:
    """
    Class for bed lines with 3 values: chrom, start and end
    """
    def __init__(self, chrom, start, end):
        self.chrom = chrom;
        self.start = start;
        self.end = end;


    def __set__(self, chrom, start, end):
        self.chrom = chrom;
        self.start = start;
        self.end = end;


    def getCoord(self):
        outstring = self.chrom + " " + str(self.start) + " " + str(self.end);
        try:
            return outstring;
        except:
            sys.stderr.write("No coord information for %s\n" % self)
            return ''


#----------------------------
#----------------------------


class BED6:
    """
    Class for bed lines with 6 values:  chrom, start, end, name, score, strand
    """
        
    def __init__(self, chrom, start, end, name, score, strand):
        self.chrom = chrom;
        self.start = start;
        self.end = end;
        self.name = name;
        self.score = score;
        self.strand = strand;

    def __set__(self, chrom, start, end, name, score, strand):
        self.chrom = chrom;
        self.start = start;
        self.end = end;
        self.name = name;
        self.score = score;
        self.strand = strand;

    def getCoord(self):
        outstring = self.chrom + "\t" + str(self.start) + "\t" + \
                    str(self.end) + "\t" + self.name + "\t" + \
                    str(self.score) + "\t" + self.strand;
        try:
            return outstring;
        except:
            sys.stderr.write("No coord information for %s\n" % self)
            return ''

#----------------------------
#----------------------------


class BED_GRAPH:
    """
    Class to deal with bed graph lines: chrom, start, end, value
    This emulates the wiggle format

    """
    def __init__(self, chrom, start, end, value=0):
        self.chrom = chrom;
        self.start = start;
        self.end = end;
        self.value = value;
    def __set__(self, chrom, start, end, value):
        self.chrom = chrom;
        self.start = start;
        self.end = end;
        self.value = value;
    def getCoord(self):
        outstring = self.chrom + " " + str(self.start) + " " + str(self.end);
        try:
            return outstring;
        except:
            sys.stderr.write("No BED coord information for %s\n" % self)
            return ''
    def getAll(self):
        outstring = self.chrom + " " + str(self.start) + " " + str(self.end) + " " + str(self.value);
        try:
            return outstring;
        except:
            sys.stderr.write("No BED all information for %s\n" % self)
            return ''


#----------------------------
#----------------------------

class BED:
    """
    Class to deal with bed files and do common operations 
    """ 
    
    def __init__(self, species="hg18", file=None, bed_type="BED3", val_threshold=0):

        """ Overload __init__ so that if a threshold is given, only
        grab bed vals that are above threshold -- This won't do
        anything different for bed3 entries because there is no value
        information in these.
              
        Reads in a bed file and builds a dictionary with chromosomes
        as keys and lists of bed elements as values. The values are
	stored as floating numbers.
        """
        
        self.bed_vals = {}

        """
        initialize a dictionary with chromosomes
        """
        for c in GenomeData.species_chroms[species]:
            self.bed_vals[c] = [];

        if(file):
            if re.match(bed_type, "BED3"):
                infile = open(file);
                for line in infile:
                    """ check to make sure not a header line """
                    if not re.match("track", line):
                        line = line.strip();
                        sline = line.split();
                        
                        if len(sline) == 3:
                            bed = BED3(sline[0], atoi(sline[1]), atoi(sline[2]));
                            self.bed_vals[sline[0]].append(bed);
                        elif len(sline) == 4:
                            if atof(sline[3]) >= val_threshold:
                                bed = BED3(sline[0], atoi(sline[1]),
                                           atoi(sline[2]));
                                self.bed_vals[sline[0]].append(bed);
                        elif len(sline) == 6:
                            if atof(sline[4]) >= val_threshold:
                                bed = BED3(sline[0], atoi(sline[1]),
                                           atoi(sline[2]));
                                self.bed_vals[sline[0]].append(bed);
                        elif len(sline) >= 3:
                            bed = BED3(sline[0], atoi(sline[1]), atoi(sline[2]));
                            self.bed_vals[sline[0]].append(bed);


            
            elif re.match(bed_type, "BED_GRAPH"):
                """ If want BED_GRAPH type """
                infile = open(file);
                for line in infile:
                    """ check to make sure not a header line """
                    if not re.match("track", line):
                        line = line.strip();
                        sline = line.split();
                        
                        if len(sline) == 3:
                            sys.stderr.write("Can't make bed_graph with only \
                            3 elements")
                            raise bedError 
                        elif len(sline) == 4:
                            if atof(sline[3]) >= val_threshold:
                                bed = BED_GRAPH(sline[0], atoi(sline[1]), atoi(sline[2]),
                                           atof(sline[3]));
                                self.bed_vals[sline[0]].append(bed);
                        elif len(sline) == 6:
                            if atof(sline[4]) >= val_threshold:
                                bed = BED_GRAPH(sline[0], atoi(sline[1]), atoi(sline[2]),
                                           atof(sline[4]));
                                self.bed_vals[sline[0]].append(bed);

            elif re.match(bed_type, "BED2"):
                """ If want BED2 type """
                infile = open(file);
                for line in infile:
                    """ check to make sure not a header line """
                    if not re.match("track", line):
                        line = line.strip();
                        sline = line.split();
                        if len(sline) == 3:
                            sys.stderr.write("Need BED6 to make BE2")
                            raise bedError
                        elif len(sline) == 4:
                            sys.stderr.write("Need BED6 to make BE2")
                            raise bedError
                        elif len(sline) == 6:
                            if atof(sline[4]) >= val_threshold:
                                if plus.match(sline[5]):
                                    bed = BED2(atoi(sline[1]), sline[5]);
                                elif minus.match(sline[5]): # The BED format is open ended, hence the real start is atoi(sline[2]) - 1
                                    bed = BED2(atoi(sline[2]) - 1, sline[5]);
                                self.bed_vals[sline[0]].append(bed);

            elif re.match(bed_type, "BED6"):
                """ If want BED6 element """
                infile = open(file);
                for line in infile:
                    """ check to make sure not a header line """
                    if not re.match("track", line):
                        line = line.strip();
                        sline = line.split();
                        if len(sline) == 3:
                            sys.stderr.write("Need BED6 to make BE6")
                            raise bedError
                        elif len(sline) == 4:
                            sys.stderr.write("Need BED6 to make BE6")
                            raise bedError
                        elif len(sline) == 6:
                            if atof(sline[4]) >= val_threshold:
                                bed = BED6(sline[0], atoi(sline[1]), atoi(sline[2]),
                                           sline[3], atof(sline[4]), sline[5]);
                                self.bed_vals[sline[0]].append(bed);

    
    def keys(self):
        """
        Return a list of the keys - duplicating the function of a dictionary
        """
        return self.bed_vals.keys()

    """
    -- BED2, BED3, BED_GRAPH all ok -- BED3 and BED_GRAPH DO NOT HAVE
    STRAND INFO SO DON'T NEEED TO WORRY ABOUT THEM HERE
    -- BED2 STRANDS WERE TAKE INTO CONSIDERATION WHEN READ IN
    -- BED6 - NEED TO USE getStarts_consider_strands
    """

    def getStarts_consider_strands(self, chrom):
        """
        Return a list of starts on a given chromosome
        """
        starts = [];
        for t in self.bed_vals[chrom]:
            if plus.match(t.strand):
                starts.append(t.start);
            elif minus.match(t.strand):
                starts.append(t.end);
        try:
            return starts;
        except:
            sys.stderr.write("Having trouble returning starts %s\n" % self)
            return ''

    def getStarts(self, chrom):
        """
        Return a list of starts on a given chromosome
        """
        starts = [];
        for t in self.bed_vals[chrom]:
            starts.append(t.start);
        try:
            return starts;
        except:
            sys.stderr.write("Having trouble returning starts %s\n" % self)
            return ''


    def getEnds(self, chrom):
        """
        Return a list of starts on a given chromosome
        """
        ends = [];
        for t in self.bed_vals[chrom]:
            ends.append(t.end);
        try:
            return ends;
        except:
            sys.stderr.write("Having trouble returning ends %s\n" % self)
            return ''

    def getChroms(self):
        """
        Return a list of all the chromosomes in order (ordered keys)
        """
        try:
            return chroms[species];
        except:
            sys.stderr.write("Can't return chromosome list\n" % self)
            return ''

   
    def getNumVals(self):
        """
        Return the number of bed vals in BED instance
        """
        num = 0;
        for c in self.bed_vals.keys():
            num += len(self.bed_vals[c]);
        return num;


    def addChrom(self, chrom, bed_list):
        if self.bed_vals.has_key(chrom) and len(self.bed_vals[chrom]) > 0:
            sys.stderr.write("chromsome %s already populated\n" % chrom)
            raise bedError
        else: self.bed_vals[chrom] = bed_list;
        

    def __del__(self):
        """
        Delete, delete;
        """
        self.bed_vals.clear()


    def __contains__(self, item):
        """
        Returns  mapping iterator
        """
        return self.bed_vals.has_key(item)

    def __iter__(self):
        """
        Returns mapping iterator
        """
        return self.bed_vals.iterkeys()

    def __len__(self):
        """
        Returns number of bed_vals
        """
        return len(self.bed_vals)

    def __delitem__(self, name):
        """
        removes a chrom if its name exists in the dictionary
        -- I guess this could possible be useful at some point
        
        """
        if self.bed_vals.has_key(name):
            del self.bed_vals[name]
        else: raise bedError
        
    def __setitem__(self, name, bedlist):
        """
        Sets a new bed value
        """
        self.bed_vals[name] = bedlist

    def __getitem__(self, name):
        """
        Returns a bed_val indexed by its name or None if no such bed_val exists
        """
        if self.bed_vals.has_key(name):
            return self.bed_vals[name]
        else: raise bedError
        


        
