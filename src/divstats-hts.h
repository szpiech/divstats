/* divstats -- htslib-backed variant reader
   Copyright (C) 2015  Zachary A Szpiech

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/
#ifndef __DIVSTATS_HTS_H__
#define __DIVSTATS_HTS_H__

#include "divstats-data.h"

//Reads VCF, bgzipped VCF and BCF in a SINGLE pass, filling the haplotype
//matrix and the map together. htslib detects the format from the file itself,
//so no flag distinguishes .vcf from .vcf.gz from .bcf.
//
//This replaces readHaplotypeDataVCF + readMapDataVCF, which between them read
//the same file four times: each counted lines, closed the file and read it
//again, and the map pass re-read CHROM and POS off lines the first pass had
//already parsed.
//
//Both output pointers are allocated here and owned by the caller; release them
//with releaseHapData and releaseMapData as before. Throws 0 on any error, to
//match the convention in divstats-data.cpp.
void readVariantDataHTS(string filename, bool HEMI,
                        HaplotypeData **hapDataOut, MapData **mapDataOut);

#endif
