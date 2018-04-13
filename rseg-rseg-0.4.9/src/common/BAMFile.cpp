/*
 * Copyright (C) 2012 University of Southern California
 *                    Andrew D Smith and Qiang Song
 * Author: Qiang Song
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include <iostream>
#include <string>

#include "GenomicRegion.hpp"
#include "BAMFile.hpp"
#include "sam.h"

using std::string;
using std::cerr;
using std::endl;

BAMFile::BAMFile()
    : filename(""), mode(""), file_handler(NULL), GOOD(false) 
{
}

BAMFile::BAMFile(const std::string &f,
                 const std::string &m)
    : filename(f), mode(m)
{
    if (mode.empty())
    {
        const string ext_name = filename.substr(filename.find_last_of('.'));
        mode = ext_name == ".bam" ? "rb" : "r";
    }
    
    if ((file_handler = samopen(filename.c_str(), mode.c_str(), NULL))
        == NULL)
    {
        cerr << "Fail to open SAM/BAM file " << filename << endl;
        exit(-1);
    }
    GOOD = true;
}

BAMFile::~BAMFile()
{
    close();
}

void
BAMFile::open(const std::string &f,
              const std::string &m)
{
    filename = f;
    mode = m;
    
    if (mode.empty())
    {
        const string ext_name = filename.substr(filename.find_last_of('.'));
        mode = ext_name == ".bam" ? "rb" : "r";
    }

    if ((file_handler = samopen(filename.c_str(), mode.c_str(), NULL))
        == NULL)
    {
        cerr << "Fail to open SAM/BAM file " << filename << endl;
        GOOD = false;
        exit(-1);
    }
    GOOD = true;
}

void
BAMFile::close()
{
    if (file_handler)
    {
        samclose(file_handler);
        file_handler = NULL;
        filename = "";
        mode = "";
        GOOD = false;
    }
}


BAMFile&
operator>>(BAMFile& in, SimpleGenomicRegion& region)
{
    bam1_t * algn_p = bam_init1();

    if (samread(in.file_handler, algn_p) >= 0)
    {
        const bam1_core_t *c = &algn_p->core;
        uint32_t *cigar = bam1_cigar(algn_p);
        size_t len = 0;
        for (size_t i = 0; i < c->n_cigar; ++i)
        {
            int op = cigar[i]&0xf;
            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
                len += cigar[i]>>4;
        }
        
        region.set_chrom(in.file_handler->header->target_name[c->tid]);
        region.set_start(c->pos);
        region.set_end(c->pos + len);
    }
    else
        in.GOOD = false;
    
    bam_destroy1(algn_p);
    
    return in;
}

BAMFile& 
operator>>(BAMFile& in, GenomicRegion& region)
{
    bam1_t * algn_p = bam_init1();

    if (samread(in.file_handler, algn_p) >= 0)
    {
        const bam1_core_t *c = &algn_p->core;
        uint32_t *cigar = bam1_cigar(algn_p);
        size_t len = 0;
        for (size_t i = 0; i < c->n_cigar; ++i)
        {
            int op = cigar[i]&0xf;
            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
                len += cigar[i]>>4;
        }

        region.set_chrom(in.file_handler->header->target_name[c->tid]);
        region.set_start(c->pos);
        region.set_end(c->pos + len);
        region.set_name(bam1_qname(algn_p));
        region.set_score(c->qual);                       
        region.set_strand(c->flag & BAM_FREVERSE ? '-' : '+'); 
    }
    else
        in.GOOD = false;
    
    bam_destroy1(algn_p);
    
    return in;
}
