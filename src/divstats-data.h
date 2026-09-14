/* selscan -- a program to calculate EHH-based scans for positive selection in genomes
   Copyright (C) 2014  Zachary A Szpiech

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

#ifndef __XP_IHH_DATA_H__
#define __XP_IHH_DATA_H__
#include <string>
#include <limits>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include <map>

using namespace std;

const double MISSING = -999;
//Tajima's D needs at least one expected segregating site. With SFS
//subsampling S is fractional, so 0 < S < 1 is reachable and D is attenuated
//by roughly sqrt(S) there; see the note in tajimaD_from_sfs. On unprojected
//data S is a whole number, so this only ever excludes S == 0.
const double MIN_SEGSITES_FOR_D = 1.0;

//S is accumulated as a sum of projected bin masses, so a site that is
//CERTAIN to stay polymorphic does not sum to exactly 1: the window's
//contribution is sum over j=1..H-1 of P(X=j), which is 1 - P(0) - P(H) only
//in exact arithmetic and lands a few ulps below it in double. Without this
//tolerance the test above rejects windows holding exactly one whole site --
//42 of them in the 200-window --winsize 2 scan of tests/data/core.vcf.gz.
//The slack is far below any attenuation worth acting on: at S = 1 - 1e-6 the
//shrinkage factor is sqrt(S) = 0.9999995.
const double SEGSITES_TOL = 1e-6;

//Undefined STATISTICS are NaN, not -999. A numeric sentinel is silently
//absorbed by anything that averages a column -- mean(), quantile(), a
//smoothing window -- turning "this window has no answer" into a plausible
//looking number roughly three orders of magnitude from any real pi or D.
//NaN propagates instead, and the token written to the file is set by
//--na-string. MISSING is still used as an integer index sentinel (s_pi,
//s_S), where only its being negative matters.
const double UNDEFINED_STAT = std::numeric_limits<double>::quiet_NaN();
const char MISSING_CHAR = '9';
const char MISSING_ALLELE = '-';
const string TPED_MISSING = "-9";
const char VCF_MISSING = '.';

//Alleles are stored 2 BITS PER SITE, four to a byte, rather than one byte
//each. There are exactly four values -- '-' missing, '0', '1', and '9', the
//fill initHaplotypeData writes before a reader overwrites it -- so they fit a
//2-bit field exactly, with no spare state to get wrong. One byte per allele
//costs nhaps*nloci: 19 MB for 400 haplotypes x 50,000 sites, but 18.6 GB for
//a 2,000-haplotype scan of 10M sites, which is the scale this matters at.
//
//The codes are assigned in ASCENDING ASCII ORDER of the characters they
//replace ('-' 45 < '0' 48 < '1' 49 < '9' 57), and the four sites in a byte
//are laid out most-significant-first. Those two properties together mean a
//byte-wise comparison of two equal-length packed rows gives the same ordering
//as a comparison of the unpacked strings. That is not decoration: hfs_window
//keys a std::map on haplotype strings and pi_k2 reads that map's order back
//out through count2hap, so the ordering is observable in --pik output.
const char ALLELE_DECODE[4] = {'-', '0', '1', '9'};

inline int alleleEncode(char c)
{
  switch (c) {
    case '-': return 0;
    case '0': return 1;
    case '1': return 2;
    default:  return 3;      //'9', the pre-read fill
  }
}

//bytes needed for one haplotype row
inline long hapRowBytes(long nloci) { return (nloci + 3) >> 2; }

//site l of row `row`, as the character it used to be stored as
inline char hapGet(const char *row, long l)
{
  return ALLELE_DECODE[(row[l >> 2] >> (6 - 2 * (l & 3))) & 3];
}

//the raw 2-bit code, for comparisons that do not need the character
inline int hapCode(const char *row, long l)
{
  return (row[l >> 2] >> (6 - 2 * (l & 3))) & 3;
}

inline void hapSet(char *row, long l, char c)
{
  int shift = 6 - 2 * (l & 3);
  row[l >> 2] = (char)((row[l >> 2] & ~(3 << shift)) | (alleleEncode(c) << shift));
}

struct HaplotypeData
{
  char **data;          //packed: hapRowBytes(nloci) bytes per haplotype
  int nhaps;
  int nloci;
};

struct MapData
{
  int *physicalPos;
  double *geneticPos;
  string *locusName;
  int nloci;
  string chr;
};

struct FreqData
{
  int *count;
  int *nmissing;
  int maxMissing;
  int nloci;
  int nhaps;
};

struct array_t
{
  double *data;
  int size;
};

struct HaplotypeFrequencySpectrum {
  map<string,int> hap2count;
  multimap<int,string> count2hap;
  int *sortedCount;
  int size;
  int numUniq;
};

struct pair_t //guess it's a quartet...
{
  int start;
  int end;
  int winStart;
  int winEnd;
};

HaplotypeFrequencySpectrum *initHaplotypeFrequencySpectrum();
void releaseHaplotypeFrequencySpectrum(HaplotypeFrequencySpectrum *data);

array_t *initArray(int size, double fill = 0);
void releaseArray(array_t* data);

//allocates the arrays and populates them with -9 or "--" depending on type
MapData *initMapData(int nloci);
void releaseMapData(MapData *data);

//allocates the arrays and populates them with MISSING
FreqData *initFreqData(int nloci);
FreqData *initFreqData(HaplotypeData* data);
void releaseFreqData(FreqData *data);

//reads in map data and also does basic checks on integrity of format
//returns a populated MapData structure if successful
//throws an exception otherwise
MapData *readMapData(string filename, int expected_loci);
MapData *readMapDataTPED(string filename, int expected_loci, int expected_haps);

//allocates the 2-d array and populated it with -9
HaplotypeData *initHaplotypeData(unsigned int nhaps, unsigned int nloci);
void releaseHapData(HaplotypeData *data);

//reads in haplotype data and also does basic checks on integrity of format
//returns a populated HaplotypeData structure if successful
//throws an exception otherwise
HaplotypeData *readHaplotypeDataTPED(string filename);

//counts the number of "fields" in a string
//where a field is defined as a contiguous set of non whitespace
//characters and fields are delimited by whitespace
int countFields(const string &str);

#endif
