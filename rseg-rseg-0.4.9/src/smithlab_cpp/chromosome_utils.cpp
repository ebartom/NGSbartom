/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "chromosome_utils.hpp"

using std::vector;
using std::string;
using std::tr1::unordered_map;

static const char *digits = "0987654321";
static const char *whitespace = " \t";

// define what to do if parts are missing: if no ':' or '-', is whole
// thing chrom?

void
parse_region_name(string region_name, 
		  string& chrom, size_t &start, size_t &end) {
  
  const size_t colon_offset = region_name.find(":");
  
  // get the chromosome
  size_t chr_offset = region_name.find_last_of(whitespace, colon_offset);
  if (chr_offset == string::npos)
    chr_offset = 0;
  else
    chr_offset += 1;
  chrom = region_name.substr(chr_offset, colon_offset - chr_offset);
  
  // get the start
  const size_t start_end = region_name.find("-", colon_offset + 1);
  const string start_string = region_name.substr(colon_offset + 1,
						 start_end - colon_offset + 1);
  start = static_cast<size_t>(atoi(start_string.c_str()));
  
  // get the end
  const size_t end_end =
    region_name.find_first_not_of(digits, start_end + 1);
  const string end_string = region_name.substr(start_end + 1,
					       end_end - start_end - 1);
  end = static_cast<size_t>(atoi(end_string.c_str()));
}


static size_t 
adjust_start_pos(const size_t orig_start, const string &chrom_name) {
  static const double LINE_WIDTH = 50.0;
  const size_t name_offset = chrom_name.length() + 2; // For the '>' and '\n';
  const size_t preceding_newlines = 
    static_cast<size_t>(std::floor(orig_start / LINE_WIDTH));
  return orig_start + preceding_newlines + name_offset;
}


static size_t 
adjust_region_size(const size_t orig_start,
		   const string &chrom_name, 
		   const size_t orig_size) {
  static const double LINE_WIDTH = 50.0;
  const size_t preceding_newlines_start = 
    static_cast<size_t>(std::floor(orig_start / LINE_WIDTH));
  const size_t preceding_newlines_end = 
    static_cast<size_t>(std::floor((orig_start + orig_size) / LINE_WIDTH));
  return (orig_size + (preceding_newlines_end - preceding_newlines_start));
}


void 
extract_regions_chrom_fasta(const string &chrom_name,
			    const string &filename, 
			    const vector<SimpleGenomicRegion> &regions,
			    vector<string> &sequences) {
  
  std::ifstream in(filename.c_str());
  for (vector<SimpleGenomicRegion>::const_iterator i(regions.begin());
       i != regions.end(); ++i) {
    
    const size_t orig_start_pos = i->get_start();
    const size_t orig_end_pos = i->get_end();
    const size_t orig_region_size = orig_end_pos - orig_start_pos;

    const size_t start_pos = adjust_start_pos(orig_start_pos, chrom_name);
    const size_t region_size = adjust_region_size(
						  orig_start_pos, chrom_name, orig_region_size);
    assert(start_pos >= 0);

    in.seekg(start_pos);
    char buffer[region_size + 1];
    buffer[region_size] = '\0';
    in.read(buffer, region_size);

    std::remove_if(buffer, buffer + region_size,
		   std::bind2nd(std::equal_to<char>(), '\n'));
    buffer[orig_region_size] = '\0';
    
    sequences.push_back(buffer);
    std::transform(sequences.back().begin(), sequences.back().end(),
		   sequences.back().begin(), std::ptr_fun(&toupper));
    assert(i->get_width() == sequences.back().length());
  }
  in.close();
}


void
extract_regions_chrom_fasta(const string &chrom_name,
			    const string &filename, 
			    const vector<GenomicRegion> &regions,
			    vector<string> &sequences) {
  
  std::ifstream in(filename.c_str());
  for (vector<GenomicRegion>::const_iterator i(regions.begin());
       i != regions.end(); ++i) {
    
    const size_t orig_start_pos = i->get_start();
    const size_t orig_end_pos = i->get_end();
    const size_t orig_region_size = orig_end_pos - orig_start_pos;

    const size_t start_pos = adjust_start_pos(orig_start_pos, chrom_name);
    const size_t region_size = adjust_region_size(
						  orig_start_pos, chrom_name, orig_region_size);
    assert(start_pos >= 0);

    in.seekg(start_pos);
    char buffer[region_size + 1];
    buffer[region_size] = '\0';
    in.read(buffer, region_size);

    std::remove_if(
		   buffer, buffer + region_size,
		   std::bind2nd(std::equal_to<char>(), '\n'));
    buffer[orig_region_size] = '\0';

    sequences.push_back(buffer);
    std::transform(
		   sequences.back().begin(), sequences.back().end(),
		   sequences.back().begin(), std::ptr_fun(&toupper));
    if (i->neg_strand())
      revcomp_inplace(sequences.back());
    assert(i->get_width() == sequences.back().length());
  }
  in.close();
}


void
extract_regions_fasta(const string &dirname,
		      const vector<GenomicRegion> &regions_in, 
		      vector<string> &sequences) {
  
  static const string FASTA_SUFFIX(".fa");
  assert(check_sorted(regions_in));
  
  vector<string> filenames;
  read_dir(dirname, filenames);

  vector<vector<GenomicRegion> > regions;
  separate_chromosomes(regions_in, regions);

  std::tr1::unordered_map<string, size_t> chrom_regions_map;
  for (size_t i = 0; i < filenames.size(); ++i)
    chrom_regions_map[strip_path(filenames[i])] = i;

  for (size_t i = 0; i < regions.size(); ++i) {
    // GET THE RIGHT FILE
    const string chrom_name(regions[i].front().get_chrom());
    const string chrom_file(chrom_name + FASTA_SUFFIX);
    std::tr1::unordered_map<string, size_t>::const_iterator f_idx =
      chrom_regions_map.find(chrom_file);
    if (f_idx == chrom_regions_map.end())
      throw SMITHLABException("chrom not found:\t" + chrom_file);
    extract_regions_chrom_fasta(
				chrom_name, filenames[f_idx->second], regions[i], sequences);
  }
}


void extract_regions_fasta(const string &dirname,
			   const vector<SimpleGenomicRegion> &regions_in, vector<string> &sequences) {

  static const string FASTA_SUFFIX(".fa");
  assert(check_sorted(regions_in));

  vector<string> filenames;
  read_dir(dirname, filenames);

  vector<vector<SimpleGenomicRegion> > regions;
  separate_chromosomes(regions_in, regions);

  std::tr1::unordered_map<string, size_t> chrom_regions_map;
  for (size_t i = 0; i < filenames.size(); ++i)
    chrom_regions_map[strip_path(filenames[i])] = i;

  for (size_t i = 0; i < regions.size(); ++i) {
    // GET THE RIGHT FILE
    const string chrom_name(regions[i].front().get_chrom());
    const string chrom_file(chrom_name + FASTA_SUFFIX);
    std::tr1::unordered_map<string, size_t>::const_iterator f_idx =
      chrom_regions_map.find(chrom_file);
    if (f_idx == chrom_regions_map.end())
      throw SMITHLABException("chrom not found:\t" + chrom_file);
    extract_regions_chrom_fasta(
				chrom_name, filenames[f_idx->second], regions[i], sequences);
  }
}


void
identify_chromosomes(const string chrom_file, const string fasta_suffix, 
		     unordered_map<string, string> &chrom_files) {
  vector<string> the_files;
  if (isdir(chrom_file.c_str())) {
    read_dir(chrom_file, fasta_suffix, the_files);
    for (size_t i = 0; i < the_files.size(); ++i)
      chrom_files[strip_path_and_suffix(the_files[i])] = the_files[i];
  }
  else chrom_files[strip_path_and_suffix(chrom_file)] = chrom_file;
}


void
identify_and_read_chromosomes(const string chrom_file, 
			      const string fasta_suffix, 
			      unordered_map<string, string> &chrom_files) {
  vector<string> the_files;
  if (isdir(chrom_file.c_str())) {
    read_dir(chrom_file, fasta_suffix, the_files);
  }
  else
    the_files.push_back(chrom_file);

  for (size_t i = 0; i < the_files.size(); ++i) {
    vector<string> names, seqs;
    read_fasta_file(the_files[i], names, seqs);
    for (size_t j = 0; j < names.size(); ++j)
      chrom_files[names[j]] = the_files[i];
  }
}
