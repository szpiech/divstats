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
#include <iostream>
#include <fstream>
#include <string>
#include "param_t.h"
#include "divstats-wintools.h"
#include "divstats-winstats.h"
#include "divstats-data.h"
#include "divstats-cli.h"

using namespace std;

int main(int argc, char *argv[])
{
  cerr << "divstats v" + VERSION + "\n";
  param_t params;
  params.setPreamble(PREAMBLE);

  params.addFlag(ARG_THREADS, DEFAULT_THREADS, "", HELP_THREADS);

  // I/O flags
  params.addFlag(ARG_FILENAME_TPED, DEFAULT_FILENAME_TPED, "", HELP_FILENAME_TPED);
  params.addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
  params.addFlag(ARG_FILENAME_POP1_VCF, DEFAULT_FILENAME_POP1_VCF, "", HELP_FILENAME_POP1_VCF);
  params.addFlag(ARG_FILENAME_MAP, DEFAULT_FILENAME_MAP, "", HELP_FILENAME_MAP);
  params.addFlag(ARG_PMAP, DEFAULT_PMAP, "", HELP_PMAP);
    
  // Window control flags
  params.addFlag(ARG_BP, DEFAULT_BP, "", HELP_BP);
  params.addFlag(ARG_SITES, DEFAULT_SITES, "", HELP_SITES);
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

  // Other flags
  params.addFlag(ARG_EHH_PART, DEFAULT_EHH_PART, "", HELP_EHH_PART);
  params.addFlag(ARG_2_SWEEPFINDER, DEFAULT_2_SWEEPFINDER, "SILENT", HELP_2_SWEEPFINDER);

  try {
    params.parseCommandLine(argc, argv);
  }
  catch (...) {
    return 1;
  }

  int numThreads = params.getIntFlag(ARG_THREADS);

  // I/O
  string tpedFilename = params.getStringFlag(ARG_FILENAME_TPED);
  bool TPED = (tpedFilename.compare(DEFAULT_FILENAME_TPED) == 0) ? false : true;
  string vcfFilename = params.getStringFlag(ARG_FILENAME_POP1_VCF);
  bool VCF = (vcfFilename.compare(DEFAULT_FILENAME_POP1_VCF) == 0) ? false : true;
  string mapFilename = params.getStringFlag(ARG_FILENAME_MAP);
  bool MAP = (mapFilename.compare(DEFAULT_FILENAME_MAP) == 0) ? false : true;
  string outfileBase = params.getStringFlag(ARG_OUTFILE);

  // Window control
  bool USE_BP = params.getBoolFlag(ARG_BP);
  bool USE_SITES = params.getBoolFlag(ARG_SITES);
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

  // Other flags
  bool SWEEPFINDER = params.getBoolFlag(ARG_2_SWEEPFINDER);
  bool EHH_PART = params.getBoolFlag(ARG_EHH_PART);
  bool PMAP = params.getBoolFlag(ARG_PMAP);

  // Check for consistency errors within flags
  bool ERROR = false;

  if (!USE_SITES && !USE_BP){
    cerr << "ERROR: Must choose to measure windows in either sites or bps.\n";
    ERROR = true;
  }

  if (USE_SITES && USE_BP){
    cerr << "ERROR: Must choose to measure windows in either sites or bps not both.\n";
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

  if(EHH_PART && (!CALC_EHH || !DO_PARTITION)) {
    cerr << "ERROR: Must specify " << ARG_EHH << " in order to use " << ARG_EHH_PART << ".\n";
    ERROR = true;
  }

  bool NEED_GMAP = (CALC_EHHK || CALC_EHH || EHH_PART) && !PMAP;

  if (numThreads <= 0) {
    cerr << "ERROR: Must specify a positive number of threads.\n";
    ERROR = true;
  }

  if (!TPED && !VCF) {
    cerr << "ERROR: Must provide a file with genetic data.\n";
    ERROR = true;
  }

  if (TPED && VCF) {
    cerr << "ERROR: Must provide a TPED or VCF not both.\n";
    ERROR = true;
  }

  if ( NEED_GMAP && !MAP ){
    cerr << "ERROR: Must provide a mapfile or set --pmap.\n";
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
  else if (PARTITIONS.size() > MAX_PARTITION) {
    cerr << "ERROR: Request for " << PARTITIONS.size() << " partitions exceeds maximum allowed (" << WINSIZE << ").\n";
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

  
  if (ERROR) {
    return 1;
  }

  string outfile = outfileBase + ".divstats.out";
  ofstream fout;
  fout.open(outfile.c_str());
  if (fout.fail()) {
    cerr << "ERROR: Failed to open " << outfile << " for writing.\n";
    return 1;
  }


  HaplotypeData *hapData;
  MapData *mapData;
  FreqData *freqData;
  if (TPED){
    hapData = readHaplotypeDataTPED(tpedFilename);
  }
  else if (VCF){
    hapData = readHaplotypeDataVCF(vcfFilename); 
  }
  if (NEED_GMAP){
    mapData = readMapData(mapFilename, hapData->nloci);
  }
  else{//load physical positions
    if(TPED) mapData = readMapDataTPED(tpedFilename, hapData->nloci, hapData->nhaps);
    else if (VCF) mapData = readMapDataVCF(vcfFilename, hapData->nloci);
  }
  freqData = initFreqData(hapData);

  if (SWEEPFINDER) {
    cout << "position\tx\tn\tfolded\n";
    for (int i = 0; i < freqData->nloci; i++) {
      cout << mapData->physicalPos[i] << "\t" << freqData->count[i] << "\t" << freqData->nhaps << "\t0\n";
    }
    return 0;
  }

  vector< pair_t* > *windows = findAllWindows(mapData, WINSIZE, WINSTEP, USE_BP);

  int numStats = (CALC_PI +
                  CALC_PIK * PIK_CHOICE.size() +
                  CALC_S +
                  CALC_TAJ_D +
                  CALC_FAY_WU_H) *
                 (DO_PARTITION * PARTITIONS.size() + 1) +
                 (EHH_PART + EHH_PART * CALC_EHHK * EHHK_CHOICES.size()) * (DO_PARTITION * PARTITIONS.size()) +
                 (CALC_EHH * EHH_WINS.size() +
                  CALC_EHHK * EHHK_CHOICES.size() * EHH_WINS.size());

  cerr << "Calculating " << numStats << " statistics in " << windows->size() << " windows.\n";

  string names = "";
  double **results = new double*[windows->size()];
  for (int i = 0; i < windows->size(); i++) results[i] = new double[numStats];

  work_order_t *order;
  pthread_t *peer = new pthread_t[numThreads];
  int prev_index = 0;
  for (int i = 0; i < numThreads; i++)
  {
    order = new work_order_t;
    order->id = i;
    order->numStats = numStats;
    order->hapData = hapData;
    order->mapData = mapData;
    order->freqData = freqData;
    //order->flog = &flog;
    //order->bar = &pbar;
    order->params = &params;
    order->results = results;
    order->windows = windows;
    order->names = &names;
    order->DO_PARTITION = DO_PARTITION;
    order->USE_BP = USE_BP;
    pthread_create(&(peer[i]),
                   NULL,
                   (void *(*)(void *))calc_stats,
                   (void *)order);
  }

  for (int i = 0; i < numThreads; i++)
  {
    pthread_join(peer[i], NULL);
  }

  delete [] peer;


  fout << "chr\tstart\tend\tnSNPs\t" << names << endl;
  for (int w = 0; w < windows->size(); w++) {
    fout << mapData->chr << "\t" 
      << windows->at(w)->winStart << "\t" 
      << windows->at(w)->winStart + WINSIZE << "\t" 
      << windows->at(w)->end - windows->at(w)->start + 1;
    for (int s = 0; s < numStats; s++) {
      fout << "\t" << results[w][s];
    }
    fout << endl;
  }

  fout.close();

  releaseHapData(hapData);
  releaseMapData(mapData);
  releaseFreqData(freqData);

  return 0;
}

