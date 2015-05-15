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
#include "divstats-winstats.h"


double pi_numerator_btw_pools(string *haps1, int length1, string *haps2, int length2, map<string, int> &hap2count) {
   double num = 0;

   for (int i = 0; i < length1; i++) {
      for (int j = 0; j < length2; j++) {
         num += hamming_dist_str(haps1[i], haps2[j]) * hap2count[haps1[i]] * hap2count[haps2[j]];
      }
   }

   return num;
}

double pi_numerator(string *haps, int length, map<string, int> &hap2count) {
   double num = 0;

   for (int i = 0; i < length; i++) {
      for (int j = i + 1; j < length; j++) {
         num += hamming_dist_str(haps[i], haps[j]) * hap2count[haps[i]] * hap2count[haps[j]];
      }
   }

   return num;
}

double pi_k2(HaplotypeFrequencySpectrum *hfs, int k) {

   pair <multimap<int, string>::iterator, multimap<int, string>::iterator> ret;
   multimap<int, string>::iterator it;

   k = (hfs->numUniq < k) ? hfs->numUniq : k;

   string *haps = new string[k];

   //counts the number of unique hap classes upto k
   int howmanyUniqHaps = 0;
   //If the kth most frequent haplotype class has > 1 haplotype associated with it, this counts how many
   int numNextClass = 0;
   //If the kth most frequent haplotype class has > 1 haplotype associated with it, this stores the ties
   string *equalFreqHaps;//length == numNextClass
   int h = 0;
   for (int i = 0; i < k; i++) {
      int key = hfs->sortedCount[i];
      numNextClass = hfs->count2hap.count(key);

      if (howmanyUniqHaps + numNextClass > k) {
         equalFreqHaps = new string[numNextClass];
         int j = 0;
         ret = hfs->count2hap.equal_range(key);
         for (it = ret.first; it != ret.second; it++) {
            equalFreqHaps[j] = it->second;
            j++;
         }
         break;
      }

      howmanyUniqHaps += numNextClass;
      ret = hfs->count2hap.equal_range(key);

      for (it = ret.first; it != ret.second; it++) {
         haps[h] = it->second;
         h++;
      }
      numNextClass = 0;
   }

   int numHapsMissing = k - howmanyUniqHaps;

   //cout << "-----\n";
   for (int i = 0; i < howmanyUniqHaps; i++) {
      //cout << "  " << haps[i] << " " << hfs->hap2count[haps[i]] << endl;
   }

   double pi = 0;
   int nhaps = 0;
   double denominator;

   if (numHapsMissing == 0) {
      for (int i = 0; i < k; i++) nhaps += hfs->hap2count[haps[i]];
      denominator = (nhaps) * (nhaps - 1) * 0.5;
      //cout << pi / denominator << endl;
      return pi_numerator(haps, k, hfs->hap2count) / denominator;
   }
   else {
      for (int i = 0; i < k; i++) {
         if (i < howmanyUniqHaps) {
            nhaps += hfs->hap2count[haps[i]];
         }
         else{
            nhaps += hfs->hap2count[equalFreqHaps[i-howmanyUniqHaps]];
         }
      }
      denominator = (nhaps) * (nhaps - 1) * 0.5;
      double pi_partial = pi_numerator(haps, howmanyUniqHaps, hfs->hap2count);
      gsl_combination * c;
      c = gsl_combination_calloc (numNextClass, numHapsMissing);

      //cout << "--next class--\n";
      string *chosenHaps = new string[numHapsMissing];
      double pi_combo = 0;
      do
      {
         //cout << "[ ";
         for (int i = 0; i < numHapsMissing; i++) {
            chosenHaps[i] = equalFreqHaps[gsl_combination_get(c, i)];
            //cout << chosenHaps[i] << "\n  ";
         }
         //cout << "] " << hfs->hap2count[chosenHaps[0]] << " : ";
         pi += pi_partial +
               pi_numerator_btw_pools(haps, howmanyUniqHaps, chosenHaps, numHapsMissing, hfs->hap2count) +
               pi_numerator(chosenHaps, numHapsMissing, hfs->hap2count);
         pi_combo = pi_partial +
                    pi_numerator_btw_pools(haps, howmanyUniqHaps, chosenHaps, numHapsMissing, hfs->hap2count) +
                    pi_numerator(chosenHaps, numHapsMissing, hfs->hap2count);
         //cout << pi_combo << " / " << denominator << " -> " << pi_combo / denominator << endl;

      } while (gsl_combination_next (c) == GSL_SUCCESS);
      pi /= nCk(numNextClass, numHapsMissing);
      gsl_combination_free (c);
      delete [] chosenHaps;
   }

   //cout << pi / denominator << endl;


   if (equalFreqHaps != NULL) {
      delete [] equalFreqHaps;
   }

   delete [] haps;

   return pi / denominator;




   //Grab the first k most frequent haplotypes
   //Right now if there are ties in the final one
   //we only take the first few upto k total haplotypes
   //In the future, will take the mean
   /*
      int i = 0;
      for (int j = 0; j < hfs->size; j++) {
         int key = hfs->sortedCount[j];
         pair <multimap<int, string>::iterator, multimap<int, string>::iterator> ret;
         ret = hfs->count2hap.equal_range(key);
         multimap<int, string>::iterator it;

         for (it = ret.first; it != ret.second; it++) {
            haps[i] = it->second;
            i++;
            if (i >= k) break;
         }
         if (it != ret.second) break;
      }

      double pi_k = 0;
      int nhaps = 0;
      for (i = 0; i < k; i++) nhaps += hfs->hap2count[haps[i]];

      double denominator = (nhaps) * (nhaps - 1) * 0.5;

      for (i = 0; i < k; i++) {
         for (int j = i + 1; j < k; j++) {
            pi_k += hamming_dist_str(haps[i], haps[j]) * hfs->hap2count[haps[i]] * hfs->hap2count[haps[j]];
         }
      }

      delete [] haps;

      return pi_k / denominator;
   */
}

double pi_k(HaplotypeFrequencySpectrum *hfs, int k) {
   k = (hfs->numUniq < k) ? hfs->numUniq : k;
   string *haps = new string[k];

   //Grab the first k most frequent haplotypes
   //Right now if there are ties in the final one
   //we only take the first few upto k total unique haplotypes
   //In the future, will take the mean
   int i = 0;
   for (int j = 0; j < hfs->size; j++) {
      int key = hfs->sortedCount[j];
      pair <multimap<int, string>::iterator, multimap<int, string>::iterator> ret;
      ret = hfs->count2hap.equal_range(key);
      multimap<int, string>::iterator it;
      for (it = ret.first; it != ret.second; it++) {
         haps[i] = it->second;
         i++;
         if (i >= k) break;
      }
      if (it != ret.second) break;
   }

   double pi_k = 0;
   int nhaps = 0;
   for (i = 0; i < k; i++) nhaps += hfs->hap2count[haps[i]];

   double denominator = (nhaps) * (nhaps - 1) * 0.5;

   for (i = 0; i < k; i++) {
      for (int j = i + 1; j < k; j++) {
         pi_k += hamming_dist_str(haps[i], haps[j]) * hfs->hap2count[haps[i]] * hfs->hap2count[haps[j]];
      }
   }

   delete [] haps;

   return pi_k / denominator;
}

HaplotypeFrequencySpectrum *hfs_window(HaplotypeData *hapData, pair_t* snpIndex) {
   HaplotypeFrequencySpectrum *hfs = initHaplotypeFrequencySpectrum();

   //Generate haplotypes and populate hap2count
   for (int hap = 0; hap < hapData->nhaps; hap++) {
      string haplotype;

      for (int site = snpIndex->start; site <= snpIndex->end; site++) {
         if (site == snpIndex->start) {
            //haplotypeList[hap] = data[hap][site];
            haplotype = hapData->data[hap][site];
         }
         else {
            //haplotypeList[hap] += data[hap][site];
            haplotype += hapData->data[hap][site];
         }
      }

      if (hfs->hap2count.count(haplotype) == 0) {
         hfs->hap2count[haplotype] = 1;
      }
      else {
         hfs->hap2count[haplotype]++;
      }
   }

   //Populate count2hap and sortedCounts
   int *sortedCount = new int[hfs->hap2count.size()]; //could contain duplicates
   hfs->numUniq = hfs->hap2count.size();
   map<string, int>::iterator it;
   int i = 0;
   for (it = hfs->hap2count.begin(); it != hfs->hap2count.end(); it++, i++) {
      sortedCount[i] = it->second;//unsorted
      hfs->count2hap.insert(pair<int, string>(it->second, it->first));
   }

   qsort(sortedCount, hfs->hap2count.size(), sizeof(int), compare);//sorted but with possible duplicates
   hfs->sortedCount = uniqInt(sortedCount, hfs->numUniq, hfs->size);//remove duplicates
   delete [] sortedCount;

   return hfs;
}

int *uniqInt(int *array, int size, int &newSize) {
   map<int, int> uniq;
   for (int i = 0; i < size; i++) {
      uniq[array[i]] = 1;
   }
   newSize = uniq.size();
   int *newArray = new int[newSize];
   int prev = array[0];
   int j = 0;
   newArray[j] = prev;
   j++;
   for (int i = 1; i < size; i++) {
      if (array[i] != prev) {
         prev = array[i];
         newArray[j] = prev;
         j++;
      }
   }
   return newArray;
}

int compare (const void *a, const void *b)
{
   return ( *(int *)b - * (int *)a );
}

double pi_window(HaplotypeData *hapData, pair_t* snpIndex) {
   //int startSnpIndex; int endSnpIndex;
   double pi = 0;
   double denominator = (hapData->nhaps) * (hapData->nhaps - 1) * 0.5;
   int length = snpIndex->end - snpIndex->start + 1;
   if (length == 0) {
      pi = 0;
   }
   else {
      for (int i = 0; i < hapData->nhaps; i++) {
         for (int j = i + 1; j < hapData->nhaps; j++) {
            pi += hamming_dist_ptr(hapData->data[i] + snpIndex->start, hapData->data[j] + snpIndex->start, length);
         }
      }
   }
   return (pi / denominator);
}

array_t *sfs_window(FreqData *freqData, pair_t* snpIndex) {

   array_t *sfs = initArray(freqData->nhaps + 1);

   for (int i = snpIndex->start; i <= snpIndex->end; i++) {
      sfs->data[freqData->count[i]]++;
   }

   return sfs;
}

double pi_from_sfs(array_t *sfs) {
   double pi = 0;
   int n = sfs->size - 1;
   double denominator = n * (n - 1) * 0.5;

   for (int i = 1; i < n; i++) {
      pi += i * (n - i) * sfs->data[i];
   }
   return pi / denominator;
}

int segsites(array_t *sfs){
   double s = 0;
   int n = sfs->size - 1;
   for (int i = 1; i < n; i++) {
      s += sfs->data[i];
   }
   return s;
}

double s_from_sfs(array_t *sfs) {
   double s = 0;
   int n = sfs->size - 1;

   for (int i = 1; i < n; i++) {
      s += sfs->data[i];
   }
   return s;
}

double a1(int n) {
   double a = 0;
   for (double i = 1; i < n; i++)
      a += 1 / i;
   return a;
}

double a2(int n) {
   double a = 0;
   for (double i = 1; i < n; i++)
      a += 1 / (i * i);
   return a;
}