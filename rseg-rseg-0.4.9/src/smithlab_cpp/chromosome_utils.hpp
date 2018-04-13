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

#ifndef CHROMOSOME_UTILS_HPP
#define CHROMOSOME_UTILS_HPP

#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iostream>
#include <ostream>
#include <limits>
#include <iterator>
#include <cassert>
#include <cmath>

#include "GenomicRegion.hpp"

void
parse_region_name(std::string region_name,
		  std::string &chrom, size_t &start, size_t &end);


// template <class T> std::string
// assemble_region_name(const T &region) {
//   return (region.get_chrom() + ":" + smithlab::toa(region.get_start()) + "-" +
//           smithlab::toa(region.get_end()));
// }


// template <class T> std::string
// assemble_region_name(const T &region, const std::string sep) {
//   return (region.get_chrom() + sep + smithlab::toa(region.get_start()) + sep +
// 	  smithlab::toa(region.get_end()));
// }


void
extract_regions_chrom_fasta(const std::string &chrom_name, 
			    const std::string &filename,
			    const std::vector<GenomicRegion> &regions, 
			    std::vector<std::string> &sequences);


void
extract_regions_chrom_fasta(const std::string &chrom_name, 
			    const std::string &filename,
			    const std::vector<SimpleGenomicRegion> &regions, 
			    std::vector<std::string> &sequences);


void
extract_regions_fasta(const std::string &dirname, 
		      const std::vector<SimpleGenomicRegion> &regions_in, 
		      std::vector<std::string> &sequences);


void
extract_regions_fasta(const std::string &dirname, 
		      const std::vector<GenomicRegion> &regions_in, 
		      std::vector<std::string> &sequences);


void
identify_chromosomes(const std::string chrom_file,
		     const std::string fasta_suffix,
		     std::tr1::unordered_map<std::string, 
					     std::string> &chrom_files);


void
identify_and_read_chromosomes(const std::string chrom_file,
			      const std::string fasta_suffix,
			      std::tr1::unordered_map<std::string,
						      std::string> &chrom_files);

#endif
