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

#include "divstats-data.h"
#include <map>
#include <vector>
#include <cstdlib>
#include "binom.h"

int compare (const void *a, const void *b);
int *uniqInt(int *array, int size, int &newSize);

double subsample_sfs(array_t *sfs, int H, int j);
void project_sfs(array_t *sub, int H, array_t *out);
array_t *sfs_window(FreqData *freqData, pair_t* snpIndex, bool SFS_SUB, bool CONST_N);
HaplotypeFrequencySpectrum *hfs_window(HaplotypeData *hapData, pair_t* snpIndex);

double pi_window(HaplotypeData *hapData, pair_t* snpIndex);
double pi_from_sfs(array_t *sfs);
double pi_k(HaplotypeFrequencySpectrum *hfs, int k);
double pi_k2(HaplotypeFrequencySpectrum *hfs, int k, pair_t *subset_snps = NULL);
double pi_numerator_btw_pools(string *haps1, int length1, string *haps2, int length2, map<string, int> &hap2count, pair_t *subset_snps = NULL);
double pi_numerator(string *haps, int length, map<string, int> &hap2count, pair_t *subset_snps = NULL);
int hamming_dist_str(string one, string two, pair_t *subset_snps = NULL);
int hamming_dist_ptr(short *one, short *two, int length);
int hamming_dist_ptr(char *one, char *two, int length);

double tajimaD_from_sfs(array_t *sfs, double pi = -9, double S = -9);
double thetaH_from_sfs(array_t *sfs);
double fayWuH_from_sfs(array_t *sfs, double pi = -9);

double ehh_from_hfs(HaplotypeFrequencySpectrum *hfs);
double ehhk_from_hfs(HaplotypeFrequencySpectrum *hfs, int k);

double segsites(array_t *sfs);
double calc_a1(int n);
double calc_a2(int n);
double calc_e1(int n, double a1);
double calc_e2(int n, double a1, double a2);

int numSitesInDataWin(pair_t* win);

//const double MISSING = -999999999;

#endif