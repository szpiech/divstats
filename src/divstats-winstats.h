/* divstats -- a program to calculate window-based diversity statistics
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

#ifndef __DIVSTATS_WINSTATS_H__
#define __DIVSTATS_WINSTATS_H__

#include "hamming_t.h"
#include "divstats-data.h"
#include <map>
#include <cstdlib>
#include "binom.h"
#include <gsl/gsl_combination.h>

int compare (const void *a, const void *b);
int *uniqInt(int *array, int size, int &newSize);
array_t *sfs_window(FreqData *freqData, pair_t* snpIndex);

HaplotypeFrequencySpectrum *hfs_window(HaplotypeData *hapData, pair_t* snpIndex);

double pi_window(HaplotypeData *hapData, pair_t* snpIndex);
double pi_from_sfs(array_t *sfs);
double pi_k(HaplotypeFrequencySpectrum *hfs, int k);
double pi_k2(HaplotypeFrequencySpectrum *hfs, int k);
double pi_numerator_btw_pools(string *haps1, int length1, string *haps2, int length2, map<string, int> &hap2count);
double pi_numerator(string *haps, int length, map<string, int> &hap2count);

double a1(int n);
double a2(int n);

#endif