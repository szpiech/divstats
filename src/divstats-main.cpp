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
#include "divstats-hts.h"
#include "divstats-winstats.h"
#include "divstats-data.h"
#include "divstats-cli.h"

#include <cmath>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[])
{
  //--version is handled before anything else, and before the banner, so that
  //`divstats --version` writes exactly one line to STDOUT and exits 0. The
  //version was previously only reachable in the stderr banner, mixed in with
  //progress output, which conda recipes and pipeline version-capture cannot
  //use. Checked here rather than through param_t because it must not require
  //the rest of a valid command line.
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == ARG_VERSION) {
      cout << VERSION << "\n";
      return 0;
    }
  }

  cerr << "divstats v" + VERSION + "\n";
  param_t params;
  params.setPreamble(PREAMBLE);

  params.addFlag(ARG_THREADS, DEFAULT_THREADS, "", HELP_THREADS);

  // I/O flags
  params.addFlag(ARG_FILENAME_TPED, DEFAULT_FILENAME_TPED, "", HELP_FILENAME_TPED);
  params.addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
  params.addFlag(ARG_FILENAME_POP1_VCF, DEFAULT_FILENAME_POP1_VCF, "", HELP_FILENAME_POP1_VCF);
  params.addFlag(ARG_FILENAME_MAP, DEFAULT_FILENAME_MAP, "", HELP_FILENAME_MAP);
  params.addFlag(ARG_HEMI, DEFAULT_HEMI, "", HELP_HEMI);
    
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
  params.addFlag(ARG_NO_SFS_SUB, DEFAULT_NO_SFS_SUB, "", HELP_NO_SFS_SUB);
  params.addFlag(ARG_CONST_N_SUB, DEFAULT_CONST_N_SUB, "", HELP_CONST_N_SUB);
  params.addFlag(ARG_2_SWEEPFINDER, DEFAULT_2_SWEEPFINDER, "", HELP_2_SWEEPFINDER);
  params.addFlag(ARG_PMAP, DEFAULT_PMAP, "", HELP_PMAP);
  params.addFlag(ARG_NA_STRING, DEFAULT_NA_STRING, "", HELP_NA_STRING);
  params.addFlag(ARG_VERSION, DEFAULT_VERSION, "", HELP_VERSION);
  params.addFlag(ARG_TARGET_N, DEFAULT_TARGET_N, "", HELP_TARGET_N);
  params.addFlag(ARG_WINDOW_N_SUB, DEFAULT_WINDOW_N_SUB, "", HELP_WINDOW_N_SUB);
  params.addListFlag(ARG_EHH_CM, DEFAULT_EHH_CM, "", HELP_EHH_CM);
  
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
  bool HEMI = params.getBoolFlag(ARG_HEMI);

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
  bool SFS_SUB = !(params.getBoolFlag(ARG_NO_SFS_SUB));

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

  //--ehh-part is validated after the --partition block below, where
  //DO_PARTITION is actually assigned. Checking it here tested DO_PARTITION
  //while it still held its initializer, so the flag could never be used.

  vector<double> EHH_CM = params.getDoubleListFlag(ARG_EHH_CM);
  bool CALC_EHH_CM = (EHH_CM[0] != 0);

  //A genetic map is needed only by --ehh-cm, which is the only thing that
  //reads geneticPos. EHH used to demand --map or --pmap and then place its
  //subwindows by physical position regardless, so the map was parsed, stored
  //and ignored -- users had to produce a recombination map to run a
  //calculation that never consulted it.
  bool NEED_GMAP = CALC_EHH_CM;

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

  if ( CALC_EHH_CM && !MAP ){
    cerr << "ERROR: --ehh-cm measures subwindows in genetic distance and needs --map.\n";
    ERROR = true;
  }

  if ( PMAP && CALC_EHH_CM ){
    cerr << "ERROR: --pmap forces physical placement, which contradicts --ehh-cm.\n";
    ERROR = true;
  }

  for (unsigned int i = 0; i < EHH_CM.size(); i++) {
    if (CALC_EHH_CM && EHH_CM[i] <= 0) {
      cerr << "ERROR: --ehh-cm widths must be > 0. Found " << EHH_CM[i] << ".\n";
      ERROR = true;
    }
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
    cerr << "ERROR: Request for " << PARTITIONS.size() << " partitions exceeds maximum allowed (" << MAX_PARTITION << ").\n";
    ERROR = true;
  }
  else {
    DO_PARTITION = true;
  }

  //Must follow the block above: DO_PARTITION is assigned there. Each
  //requirement is reported separately so the message names the flag that is
  //actually missing.
  if (EHH_PART && !CALC_EHH) {
    cerr << "ERROR: Must specify " << ARG_EHH << " in order to use " << ARG_EHH_PART << ".\n";
    ERROR = true;
  }
  if (EHH_PART && !DO_PARTITION) {
    cerr << "ERROR: Must specify " << ARG_PARTITION << " in order to use " << ARG_EHH_PART << ".\n";
    ERROR = true;
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
  //The readers signal failure with `throw 0`, and only parseCommandLine was
  //wrapped -- so a malformed input aborted with SIGABRT and "terminating due
  //to uncaught exception" after its own diagnostic. Catch here and exit 1.
  //VCF/BCF is read in one pass, filling genotypes and positions together;
  //TPED still uses the two-function path.
  MapData *vcfMap = NULL;
  try {
    if (TPED){
      hapData = readHaplotypeDataTPED(tpedFilename);
    }
    else if (VCF){
      readVariantDataHTS(vcfFilename, HEMI, &hapData, &vcfMap);
    }
    if (NEED_GMAP){
      //an explicit --map overrides the positions carried in the variant file
      if (vcfMap != NULL) releaseMapData(vcfMap);
      mapData = readMapData(mapFilename, hapData->nloci);
    }
    else{//load physical positions
      if(TPED) mapData = readMapDataTPED(tpedFilename, hapData->nloci, hapData->nhaps);
      else if (VCF) mapData = vcfMap;
    }
  }
  catch (...) {
    return 1;
  }
  freqData = initFreqData(hapData);

  //---- resolve the sample size the whole run will project to -------------
  //Before 2.0.0 each window was projected to its own minimum observed sample
  //size, so n varied from window to window with local missingness and pi, S,
  //D and H were not comparable between windows -- with n absent from the
  //output, nothing downstream could even detect it. One value is now chosen
  //for the entire run.
  //
  //The default is the largest n every site can reach, nhaps - maxMissing,
  //which discards no site. That value is set by the single worst-covered site
  //in the file, so on real data it can be far below the typical site's n; the
  //report below shows what raising it with --target-n would cost in sites, so
  //the trade is visible rather than buried.
  int TARGET_N = params.getIntFlag(ARG_TARGET_N);
  bool WINDOW_N = params.getBoolFlag(ARG_WINDOW_N_SUB);
  bool TARGET_N_EXPLICIT = (TARGET_N > 0);
  int globalMaxN = freqData->nhaps - freqData->maxMissing;

  if (WINDOW_N && TARGET_N_EXPLICIT) {
    cerr << "ERROR: --window-n-sub and --target-n choose the sample size differently.\n";
    return 1;
  }
  if (WINDOW_N) {
    TARGET_N = 0;   //per-window minimum
  }
  else if (!TARGET_N_EXPLICIT) {
    TARGET_N = globalMaxN;
  }

  if (SFS_SUB && !WINDOW_N) {
    if (TARGET_N > freqData->nhaps) {
      cerr << "ERROR: --target-n " << TARGET_N << " exceeds the " << freqData->nhaps
           << " haplotypes in the data.\n";
      return 1;
    }
    long usable = 0;
    for (int i = 0; i < freqData->nloci; i++) {
      if (freqData->nhaps - freqData->nmissing[i] >= TARGET_N) usable++;
    }
    cerr << "Projecting every window to n = " << TARGET_N << " haplotypes ("
         << usable << " of " << freqData->nloci << " sites usable"
         << (TARGET_N_EXPLICIT ? "" : ", the largest n that excludes no site") << ").\n";

    if (!TARGET_N_EXPLICIT && freqData->nloci > 0) {
      //What would a slightly higher n cost? Sort the per-site sample sizes and
      //read off the n reachable by 99.9%, 99% and 95% of sites. If one bad
      //site is holding the whole run down this makes it obvious immediately.
      vector<int> perSite(freqData->nloci);
      for (int i = 0; i < freqData->nloci; i++) perSite[i] = freqData->nhaps - freqData->nmissing[i];
      sort(perSite.begin(), perSite.end());
      const double keep[3] = {0.999, 0.99, 0.95};
      bool worth = false;
      int lastN = -1;
      string line = "  raising it would cost sites:";
      for (int k = 0; k < 3; k++) {
        int idx = (int)((1.0 - keep[k]) * freqData->nloci);
        if (idx >= freqData->nloci) idx = freqData->nloci - 1;
        int n_k = perSite[idx];
        if (n_k > TARGET_N && n_k != lastN) {
          worth = true;
          lastN = n_k;
          char buf[96];
          snprintf(buf, sizeof(buf), "  n=%d keeps %.1f%% of sites;", n_k, 100.0 * keep[k]);
          line += buf;
        }
      }
      if (worth) cerr << line << " set with --target-n.\n";
    }
  }
  else if (WINDOW_N) {
    cerr << "Projecting each window to its own minimum sample size; "
         << "windows are NOT comparable.\n";
  }

  if (SWEEPFINDER) {
    //n is the number of haplotypes ACTUALLY OBSERVED at each site. This used
    //to write freqData->nhaps -- the full sample size -- at every site,
    //regardless of missing genotypes, so a site called in half the cohort was
    //declared to SweepFinder2 as a full-depth observation. The derived count x
    //was correct, so the reported frequency x/n was biased downwards by
    //exactly the local missingness, and the likelihood surface with it.
    //
    //The file also went to stdout, which the banner and progress lines do not,
    //so it could only be captured by redirecting a stream the program also
    //uses for nothing else in this mode. It now goes to <out>.sweepfinder.out
    //alongside the other output.
    string sfFilename = outfileBase + ".sweepfinder.out";
    ofstream sfout(sfFilename.c_str());
    if (sfout.fail()) {
      cerr << "ERROR: Failed to open " << sfFilename << " for writing.\n";
      return 1;
    }
    cerr << "Writing SweepFinder2 input to " << sfFilename << "\n";

    //folded is 0, i.e. the spectrum is unfolded and the ALT allele is assumed
    //to be the derived one. divstats has no outgroup information, so this is
    //an assumption about the input, not something it can verify.
    sfout << "position\tx\tn\tfolded\n";
    long nskipped = 0;
    for (int i = 0; i < freqData->nloci; i++) {
      int n = freqData->nhaps - freqData->nmissing[i];
      if (n <= 0) { nskipped++; continue; }   //no observed haplotype here
      sfout << mapData->physicalPos[i] << "\t" << freqData->count[i]
            << "\t" << n << "\t0\n";
    }
    sfout.close();
    if (nskipped > 0) {
      cerr << "Skipped " << nskipped << " site(s) with no called genotype.\n";
    }
    return 0;
  }

  vector< pair_t* > *windows = findAllWindows(mapData, WINSIZE, WINSTEP, USE_BP);

  //Column names and the column count both come from the command line, via one
  //function that mirrors the order calc_stats fills results[][]. numStats used
  //to be computed here by arithmetic replicating calc_stats' branch structure,
  //and the names were accumulated inside a worker thread while it processed
  //window 0 -- so a run with zero windows wrote a header with no statistic
  //columns, and the two copies of the same knowledge could drift apart.
  vector<string> colNames = buildColumnNames(&params, DO_PARTITION);
  int numStats = (int)colNames.size();

  string NA_STRING = params.getStringFlag(ARG_NA_STRING);
  cerr << "Calculating " << numStats << " statistics in " << windows->size() << " windows.\n";

  //per-window sample size actually used, and the number of sites that fed the
  //spectrum -- the two things U6 says the output must carry if windows are to
  //be interpreted at all
  int *nhapsUsed  = new int[windows->size()];
  int *nSitesUsed = new int[windows->size()];
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
    order->TARGET_N = TARGET_N;
    order->nhapsUsed = nhapsUsed;
    order->nSitesUsed = nSitesUsed;
    order->hapData = hapData;
    order->mapData = mapData;
    order->freqData = freqData;
    //order->flog = &flog;
    //order->bar = &pbar;
    order->params = &params;
    order->results = results;
    order->windows = windows;
    order->DO_PARTITION = DO_PARTITION;
    order->USE_BP = USE_BP;
    order->SFS_SUB = SFS_SUB;
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
  //nhapsUsed/nSitesUsed are freed after the write loop below


  //Data rows are tab-separated, but the statistic names used to be joined with
  //spaces and appended after a single tab -- so the header had 6 tab-delimited
  //fields where the rows had 9, and both read.table(header=TRUE) and
  //pandas.read_csv(sep='\t') mis-aligned. There was a trailing space too.
  //nhaps is the sample size the window's statistics refer to. Without it,
  //output from different runs -- or from before 2.0.0, where it varied by
  //window -- cannot be told apart or pooled safely. nSNPsUsed appears only
  //when --target-n was given explicitly, because that is the only setting
  //that can exclude sites: the default target is reachable by every site.
  fout << "chr\tstart\tend\tnbps\tnSNPs\tnhaps";
  if (TARGET_N_EXPLICIT) fout << "\tnSNPsUsed";
  for (unsigned int i = 0; i < colNames.size(); i++) fout << "\t" << colNames[i];
  fout << "\n";
  for (int w = 0; w < windows->size(); w++) {
    fout << mapData->chr << "\t" 
      << windows->at(w)->winStart << "\t" 
      << windows->at(w)->winEnd << "\t"
      << windows->at(w)->winEnd - windows->at(w)->winStart + 1 << "\t"
      << windows->at(w)->end - windows->at(w)->start + 1
      << "\t" << nhapsUsed[w];
    if (TARGET_N_EXPLICIT) fout << "\t" << nSitesUsed[w];
    for (int s = 0; s < numStats; s++) {
      //An undefined statistic is NaN internally; what reaches the file is the
      //--na-string token. The default, "nan", is what iostream would print
      //anyway, so the branch matters only when the user asks for something
      //else (NA for R, an empty field, or -999 to restore 1.x output).
      if (std::isnan(results[w][s])) fout << "\t" << NA_STRING;
      else fout << "\t" << results[w][s];
    }
    fout << endl;
  }

  fout.close();
  delete [] nhapsUsed;
  delete [] nSitesUsed;

  releaseHapData(hapData);
  releaseMapData(mapData);
  releaseFreqData(freqData);

  return 0;
}

