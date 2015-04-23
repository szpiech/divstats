/* divstats -- a program to calculate window-based diversity statistics
   Copyright (C) 2015  Zachary A Szpiech

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/
#include <iostream>
#include <fstream>
#include <string>
#include "param_t.h"
#include "hamming_t.h"
#include "divstats-data.h"

using namespace std;

const string PREAMBLE = "";

// I/O flags
const string ARG_FILENAME_TPED = "--tped";
const string DEFAULT_FILENAME_TPED = "__hapfile1";
const string HELP_FILENAME_TPED = "A TPED file containing haplotype and map data.\n\
\tVariants should be coded 0/1";

const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "The basename for all output files.";

// Window control flags
const string ARG_WINSIZE = "--winsize";
const int DEFAULT_WINSIZE = 100000;
const string HELP_WINSIZE = "The window size within which to calculate diversity statistics.";

const string ARG_WINSTEP = "--winstep";
const int DEFAULT_WINSTEP = 100000;
const string HELP_WINSTEP = "The sliding window step size.";

const string ARG_PARTITION = "--partition";
const int DEFAULT_PARTITION = 0;
const string HELP_PARTITION = "Partition the sliding window into non-overlapping sub windows\n\
\tof varying sizes within which all statistics (except EHH-based ones) are calculated separately.\n\
\te.g. For a sliding window of 100kb, --partition 25000 50000 25000 would instruct\n\
\tdivstats to calculate statistics separately within a central 50kb and within the\n\
\ttwo flanking 25kb regions for each 100kb window. Partitions must add up to --winsize.\n\
\tSet to 0 to simply calculate within the entire window.";

// Statistics flags
const string ARG_PI = "--pi";
const bool DEFAULT_PI = false;
const string HELP_PI = "Set this flag to calculate mean pairwise sequence difference.";

const string ARG_PIK = "--pik";
const int DEFAULT_PIK = 0;
const string HELP_PIK = "Set this flag to calculate mean pairwise sequence difference amongst \n\
\tthe k most frequent haplotypes. You can choose more than one, e.g. --pik 2 3 4 will\n\
\tcalculate pi amongst the top 2, 3, and 4 most frequent haplotypes.  If set to 0, does\n\
\tnot calculate.";

const string ARG_SEGSITES = "--s";
const bool DEFAULT_SEGSITES = false;
const string HELP_SEGSITES = "Set this flag to calculate the number of segregating sites.";

const string ARG_EHH = "--ehh";
const int DEFAULT_EHH = 0;
const string HELP_EHH = "A list of window sizes within which to calculate EHH. These will be\n\
\tcentered on the middle of the current window and may not be larger than --winsize.\n\
\tSet to 0 to simply calculate within the entire window.";

const string ARG_EHHK = "--ehhk";
const int DEFAULT_EHHK = 0;
const string HELP_EHHK = "Calculates EHH, after collapsing\n\
\tthe k most frequent haplotypes into a single identity class using the windows\n\
\tdefined by --ehh. This flag requires --ehh to be set.\n\
\tIf set to 0 does not calculate.";

const string ARG_TAJ_D = "--d";
const bool DEFAULT_TAJ_D = false;
const string HELP_TAJ_D = "Set this flag to calculate Tajima's D.";

const string ARG_FAY_WU_H = "--h";
const bool DEFAULT_FAY_WU_H = false;
const string HELP_FAY_WU_H = "Set this flag to calculate Fay and Wu's H.";

pair_t* findInclusiveSNPIndicies(int startSnpIndex, int currWinStart, int WINSIZE, MapData* mapData);

double pi_window(HaplotypeData *hapData, pair_t* snpIndex);

int main(int argc, char *argv[])
{
  param_t params;
  params.setPreamble(PREAMBLE);

  // I/O flags
  params.addFlag(ARG_FILENAME_TPED, DEFAULT_FILENAME_TPED, "", HELP_FILENAME_TPED);
  params.addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);

  // Window control flags
  params.addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "", HELP_WINSIZE);
  params.addFlag(ARG_WINSTEP, DEFAULT_WINSTEP, "", HELP_WINSTEP);
  params.addListFlag(ARG_PARTITION, DEFAULT_PARTITION, "", HELP_PARTITION);

  // Statistics flags
  params.addFlag(ARG_PI, DEFAULT_PI, "", HELP_PI);
  params.addListFlag(ARG_PIK, DEFAULT_PIK, "", HELP_PIK);
  params.addFlag(ARG_SEGSITES, DEFAULT_SEGSITES, "", HELP_SEGSITES);
  params.addListFlag(ARG_EHH, DEFAULT_EHH, "", HELP_EHH);
  params.addListFlag(ARG_EHHK, DEFAULT_EHHK, "", HELP_EHHK);
  params.addFlag(ARG_TAJ_D, DEFAULT_TAJ_D, "", HELP_TAJ_D);
  params.addFlag(ARG_FAY_WU_H, DEFAULT_FAY_WU_H, "", HELP_FAY_WU_H);

  try {
    params.parseCommandLine(argc, argv);
  }
  catch (...) {
    return 1;
  }

  // I/O
  string tpedFilename = params.getStringFlag(ARG_FILENAME_TPED);
  string outfileBase = params.getStringFlag(ARG_OUTFILE);

  // Window control
  int WINSIZE = params.getIntFlag(ARG_WINSIZE);
  int WINSTEP = params.getIntFlag(ARG_WINSTEP);
  vector<int> PARTITIONS = params.getIntListFlag(ARG_PARTITION);
  bool DO_PARTITION = false;

  // Statistics
  bool CALC_PI = params.getBoolFlag(ARG_PI);
  vector<int> PIK_CHOICE = params.getIntListFlag(ARG_PIK);
  bool CALC_PIK = false;
  bool CALC_S = params.getBoolFlag(ARG_SEGSITES);
  vector<int> EHH_WINS = params.getIntListFlag(ARG_EHH);
  bool CALC_EHH = false;
  vector<int> EHHK_CHOICES = params.getIntListFlag(ARG_EHHK);
  bool CALC_EHHK = false;
  bool CALC_TAJ_D = params.getBoolFlag(ARG_TAJ_D);
  bool CALC_FAY_WU_H = params.getBoolFlag(ARG_FAY_WU_H);

  // Check for consistency errors within flags
  bool ERROR = false;

  if (tpedFilename.compare(DEFAULT_FILENAME_TPED) == 0) {
    cerr << "ERROR: Must provide a TPED file.\n";
    ERROR = true;
  }

  if (WINSIZE < 1) {
    cerr << "ERROR: Window size needs to be greater than 0.\n";
    ERROR = true;
  }

  if (WINSTEP < 1) {
    cerr << "ERROR: Window step size needs to be greater than 0.\n";
    ERROR = true;
  }

  int partitionTotalSize = 0;
  for (int i = 0; i < PARTITIONS.size(); i++) {
    partitionTotalSize += PARTITIONS[i];
    if (PARTITIONS.size() > 1 && PARTITIONS[i] <= 0) {
      cerr << "ERROR: Partitions must be > 0. Found partition " << i + 1 << " equal to " << PARTITIONS[i] << ".\n";
      ERROR = true;
    }
  }

  if (partitionTotalSize <= 0 && PARTITIONS.size() == 1) {
    DO_PARTITION = false;
  }
  else if (partitionTotalSize < WINSIZE || partitionTotalSize > WINSIZE) {
    cerr << "ERROR: Window partitions sum to " << partitionTotalSize << " but must sum to " << WINSIZE << " instead.\n";
    ERROR = true;
  }
  else {
    DO_PARTITION = true;
  }

  for (int i = 0; i < PIK_CHOICE.size(); i++) {
    if (PIK_CHOICE.size() > 1 && PIK_CHOICE[i] <= 0) {
      cerr << "ERROR: When chosing one or more k for the k most frequent haplotypes,\n\
      \tk must be > 0. Found k number " << i + 1 << " equal to " << PIK_CHOICE[i] << ".\n";
      ERROR = true;
    }
  }

  if (PIK_CHOICE.size() == 1 && PIK_CHOICE[0] <= 0) {
    CALC_PIK = false;
  }
  else {
    CALC_PIK = true;
  }

  for (int i = 0; i < EHH_WINS.size(); i++) {
    if (EHH_WINS.size() > 1 && (EHH_WINS[i] <= 0 || EHH_WINS[i] > WINSIZE)) {
      cerr << "ERROR: When chosing EHH windows, each must be > 0 and <= " << WINSIZE << ".\n\
      \tFound window " << i + 1 << " equal to " << EHH_WINS[i] << ".\n";
      ERROR = true;
    }
    else if (EHH_WINS.size() == 1 && EHH_WINS[i] > WINSIZE) {
      cerr << "ERROR: When chosing EHH windows, each must be > 0 and <= " << WINSIZE << ".\n\
      \tFound window " << i + 1 << " equal to " << EHH_WINS[i] << ".\n";
      ERROR = true;
    }
  }

  if (EHH_WINS.size() == 1 && EHH_WINS[0] <= 0) {
    CALC_EHH = false;
  }
  else {
    CALC_EHH = true;
  }

  for (int i = 0; i < EHHK_CHOICES.size(); i++) {
    if (EHHK_CHOICES.size() > 1 && EHHK_CHOICES[i] <= 0) {
      cerr << "ERROR: When chosing one or more k for the k most frequent haplotypes,\n\
      \tk must be > 0. Found k number " << i + 1 << " equal to " << EHHK_CHOICES[i] << ".\n";
      ERROR = true;
    }
  }

  if (EHHK_CHOICES.size() == 1 && EHHK_CHOICES[0] <= 0) {
    CALC_EHHK = false;
  }
  else {
    CALC_EHHK = true;
  }

  if (CALC_EHHK && !CALC_EHH) {
    cerr << "ERROR: Must specify " << ARG_EHH << " in order to use " << ARG_EHHK << ".\n";
    ERROR = true;
  }

  if (ERROR) {
    return 1;
  }


  HaplotypeData *hapData;
  MapData *mapData;
  FreqData *freqData;

  hapData = readHaplotypeDataTPED(tpedFilename);
  mapData = readMapDataTPED(tpedFilename, hapData->nloci, hapData->nhaps);
  freqData = initFreqData(hapData);



  int currWinStart = 1;//mapData->physicalPos[0];
  int currWinEnd = currWinStart + WINSIZE - 1;
  int numSnps = mapData->nloci;
  int endOfData = mapData->physicalPos[numSnps - 1];
  pair_t *snpIndex = new pair_t;
  snpIndex->start = 0;
  snpIndex->end = -1;
  int numInWindow;

  for (currWinStart; currWinStart < endOfData; currWinStart += WINSTEP, currWinEnd += WINSTEP) {
    vector< pair_t* > *windows = new vector< pair_t* >;

    //Find SNP index boundaries for the whole window
    pair_t *snps = findInclusiveSNPIndicies(snpIndex->start, currWinStart, WINSIZE, mapData);
    windows->push_back(snps);

    snpIndex->start = snps->start;
    snpIndex->end = snps->end;

    //Find SNP index boundaries for partitions
    if (DO_PARTITION) {
      pair_t *partitionSnpIndex = new pair_t;
      partitionSnpIndex->start = snpIndex->start;
      partitionSnpIndex->end = snpIndex->start - 1;
      int partitionCurrWinStart = currWinStart;
      for (int i = 0; i < PARTITIONS.size(); i++) {
        pair_t *partition_snps = findInclusiveSNPIndicies(partitionSnpIndex->start, partitionCurrWinStart, PARTITIONS[i], mapData);
        windows->push_back(partition_snps);
        partitionSnpIndex->start = partition_snps->start;
        partitionSnpIndex->end = partition_snps->end;
        partitionCurrWinStart += PARTITIONS[i];
      }
      delete partitionSnpIndex;
    }

    cout << currWinStart << " " << currWinEnd;
    //Cycle over all windows and partitions
    for (int i = 0; i < windows->size(); i++) {
      snps = windows->at(i);

      cout << " " << snps->end - snps->start + 1 << " " << pi_window(hapData, snps);

      snps = NULL;
      delete windows->at(i);
    }
    cout << "\n";
    delete windows;
  }

  delete snpIndex;
  return 0;
}

pair_t* findInclusiveSNPIndicies(int startSnpIndex, int currWinStart, int WINSIZE, MapData* mapData) {

  int currWinEnd = currWinStart + WINSIZE - 1;
  int endSnpIndex = startSnpIndex;
  int numSnps = mapData->nloci;

  pair_t* snps = new pair_t;

  if (mapData->physicalPos[numSnps - 1] < currWinStart) {
    snps->start = numSnps;
    snps->end = numSnps - 1;
    return snps;
  }

  while (mapData->physicalPos[startSnpIndex] < currWinStart) {
    startSnpIndex++;
  }
  while (mapData->physicalPos[endSnpIndex] < currWinEnd) {
    endSnpIndex++;
  }
  endSnpIndex--;
  endSnpIndex = (endSnpIndex >= numSnps) ? numSnps - 1 : endSnpIndex;

  snps->start = startSnpIndex;
  snps->end = endSnpIndex;
  return snps;
}

double pi_window(HaplotypeData *hapData, pair_t* snpIndex) {
  //int startSnpIndex; int endSnpIndex;
  double pi = 0;
  double denominator = (hapData->nhaps) * (hapData->nhaps - 1) * 0.5;
  int length = snpIndex->end - snpIndex->start + 1;
  if (length == 0)
  {
    pi = 0;
  }
  else
  {
    for (int i = 0; i < hapData->nhaps; i++)
    {
      for (int j = i + 1; j < hapData->nhaps; j++)
      {
        pi += hamming_dist_ptr(hapData->data[i] + snpIndex->start, hapData->data[j] + snpIndex->start, length);
      }
    }
  }
  return (pi / denominator);
}

