#include "divstats-wintools.h"
#include <cstdio>


vector< pair_t* > *findAllWindows(MapData *mapData, int WINSIZE, int WINSTEP, bool USE_BP) {
	vector< pair_t* > *windows = new vector< pair_t* >;
	int numSnps = mapData->nloci;
	if (USE_BP){
		int endOfData = mapData->physicalPos[numSnps - 1];	
		int snpIndexStart = 0;

		for (int currWinStart = 0; currWinStart < endOfData; currWinStart += WINSTEP/*, currWinEnd += WINSTEP*/) {	
			if(currWinStart+WINSIZE-1 >= endOfData) break;
			//Find SNP index boundaries for the whole window
			pair_t *snps = findInclusiveSNPIndicies(snpIndexStart, currWinStart, WINSIZE, mapData);
			windows->push_back(snps);
			snpIndexStart = snps->start;
		}
	}
	else{//USE_SITES
		for (int i = 0; i < numSnps; i += WINSTEP){
			if (i+WINSIZE-1 >= numSnps) break;
			pair_t* snps = new pair_t;
			snps->start = i;
			snps->end = i+WINSIZE-1;
			snps->winStart = mapData->physicalPos[i];
			snps->winEnd = mapData->physicalPos[i+WINSIZE-1];
			windows->push_back(snps);
		}
	}
	return windows;
}

void releaseAllWindows(vector< pair_t* > *windows) {
	for (int i = 0; i < windows->size(); i++) delete windows->at(i);
	delete windows;
	return;
}

//The column names are a pure function of the command line, but they used to be
//assembled as a side effect of thread 0 processing window 0
//(`if (i == 0) (*names) += "pi ";`). Two consequences: a run that produced zero
//windows emitted a header with no statistic columns at all, and any change to
//how windows are scheduled would silently desynchronise the header from the
//data. numStats was computed separately in main by arithmetic that mirrored the
//branch structure here by hand, which is a second copy of the same knowledge.
//
//This returns the names in the order calc_stats fills results[][], so main can
//write the header before any thread starts and take numStats from the count.
//It must stay in lockstep with the two blocks below; the runtime check in main
//compares the count against what calc_stats actually wrote.
vector<string> buildColumnNames(param_t *params, bool DO_PARTITION) {
	vector<string> names;
	vector<int> PIK_CHOICE    = params->getIntListFlag(ARG_PIK);
	vector<int> EHH_WINS      = params->getIntListFlag(ARG_EHH);
	vector<int> EHHK_CHOICES  = params->getIntListFlag(ARG_EHHK);
	vector<int> PARTITIONS    = params->getIntListFlag(ARG_PARTITION);
	vector<double> EHH_CM     = params->getDoubleListFlag(ARG_EHH_CM);
	bool EHH_PART = params->getBoolFlag(ARG_EHH_PART);

	//whole-window statistics, in STATS order
	for (int j = 0; j < NOPTS; j++) {
		if (STATS[j].compare(ARG_PI) == 0 && params->getBoolFlag(ARG_PI)) {
			names.push_back("pi");
		}
		else if (STATS[j].compare(ARG_PIK) == 0 && PIK_CHOICE[0] != 0) {
			for (unsigned int k = 0; k < PIK_CHOICE.size(); k++)
				names.push_back("pi" + int2str(PIK_CHOICE[k]));
		}
		else if (STATS[j].compare(ARG_SEGSITES) == 0 && params->getBoolFlag(ARG_SEGSITES)) {
			names.push_back("S");
		}
		else if (STATS[j].compare(ARG_EHH) == 0 && EHH_WINS[0] != 0) {
			for (unsigned int w = 0; w < EHH_WINS.size(); w++)
				names.push_back("ehh_" + int2str(EHH_WINS[w]));
		}
		else if (STATS[j].compare(ARG_EHHK) == 0 && EHHK_CHOICES[0] != 0) {
			for (unsigned int w = 0; w < EHH_WINS.size(); w++)
				for (unsigned int k = 0; k < EHHK_CHOICES.size(); k++)
					names.push_back("ehh" + int2str(EHHK_CHOICES[k]) + "_" + int2str(EHH_WINS[w]));
		}
		else if (STATS[j].compare(ARG_EHH_CM) == 0 && EHH_CM[0] != 0) {
			for (unsigned int w = 0; w < EHH_CM.size(); w++)
				names.push_back("ehhcm_" + dbl2str(EHH_CM[w]));
		}
		else if (STATS[j].compare(ARG_TAJ_D) == 0 && params->getBoolFlag(ARG_TAJ_D)) {
			names.push_back("D");
		}
		else if (STATS[j].compare(ARG_FAY_WU_H) == 0 && params->getBoolFlag(ARG_FAY_WU_H)) {
			names.push_back("H");
		}
	}

	//per-partition statistics, labelled A, B, C ... in partition order
	if (DO_PARTITION) {
		for (unsigned int p = 0; p < PARTITIONS.size(); p++) {
			string partStr(1, (char)('A' + p));
			for (int j = 0; j < NOPTS; j++) {
				if (STATS[j].compare(ARG_PI) == 0 && params->getBoolFlag(ARG_PI)) {
					names.push_back("pi_" + partStr);
				}
				else if (STATS[j].compare(ARG_PIK) == 0 && PIK_CHOICE[0] != 0) {
					for (unsigned int k = 0; k < PIK_CHOICE.size(); k++)
						names.push_back("pi" + int2str(PIK_CHOICE[k]) + "_" + partStr);
				}
				else if (STATS[j].compare(ARG_SEGSITES) == 0 && params->getBoolFlag(ARG_SEGSITES)) {
					names.push_back("S_" + partStr);
				}
				else if (STATS[j].compare(ARG_TAJ_D) == 0 && params->getBoolFlag(ARG_TAJ_D)) {
					names.push_back("D_" + partStr);
				}
				else if (STATS[j].compare(ARG_FAY_WU_H) == 0 && params->getBoolFlag(ARG_FAY_WU_H)) {
					names.push_back("H_" + partStr);
				}
				else if (STATS[j].compare(ARG_EHH) == 0 && EHH_WINS[0] != 0 && EHH_PART) {
					names.push_back("ehh_" + partStr);
				}
				else if (STATS[j].compare(ARG_EHHK) == 0 && EHHK_CHOICES[0] != 0 && EHH_PART) {
					for (unsigned int k = 0; k < EHHK_CHOICES.size(); k++)
						names.push_back("ehh" + int2str(EHHK_CHOICES[k]) + "_" + partStr);
				}
			}
		}
	}

	return names;
}

void calc_stats(void *order) {
	work_order_t *p = (work_order_t *)order;
	HaplotypeData *hapData = p->hapData;
	MapData *mapData = p->mapData;
	FreqData *freqData = p->freqData;
	param_t *params = p->params;
	vector< pair_t* > *windows = p->windows;
	int WINSIZE = params->getIntFlag(ARG_WINSIZE);
	double **results = p->results;
	int id = p->id;
	int numStats = p->numStats;
	vector<int> PIK_CHOICE = params->getIntListFlag(ARG_PIK);
	vector<int> EHH_WINS = params->getIntListFlag(ARG_EHH);
	vector<int> EHHK_CHOICES = params->getIntListFlag(ARG_EHHK);
	vector<int> PARTITIONS = params->getIntListFlag(ARG_PARTITION);
		vector<double> EHH_CM = params->getDoubleListFlag(ARG_EHH_CM);
	bool DO_PARTITION = p->DO_PARTITION;
	bool USE_BP = p->USE_BP;
	bool SFS_SUB = p->SFS_SUB;
	bool CONST_N = params->getBoolFlag(ARG_CONST_N_SUB);

	int numThreads = params->getIntFlag(ARG_THREADS);
	//These must be initialized: each is assigned only inside a conditional
	//branch below but released unconditionally at the end of every window
	//iteration, so leaving them indeterminate frees a garbage pointer.
	array_t *sfs = NULL, *partition_sfs = NULL;
	HaplotypeFrequencySpectrum *hfs = NULL, *partition_hfs = NULL, *pik_hfs = NULL;
	pair_t *snps = NULL, *partition_snps = NULL;


	bool NEED_SFS = false;
	//Do we need to calculate the SFS for every window?
	for (int j = 0; j < NOPTS; j++){
		if (STATS[j].compare(ARG_PI) == 0 && params->getBoolFlag(ARG_PI)) NEED_SFS = true;
		else if (STATS[j].compare(ARG_SEGSITES) == 0 && params->getBoolFlag(ARG_SEGSITES)) NEED_SFS = true;
		else if (STATS[j].compare(ARG_TAJ_D) == 0 && params->getBoolFlag(ARG_TAJ_D)) NEED_SFS = true;
		else if (STATS[j].compare(ARG_FAY_WU_H) == 0 && params->getBoolFlag(ARG_FAY_WU_H)) NEED_SFS = true;
	}

	//Cycle over all windows and calculate stats
	for (int i = id; i < windows->size(); i += numThreads) {
		snps = windows->at(i);

		//Reset per-window: whatever the previous iteration allocated has
		//already been released, so these must not be carried over.
		sfs = NULL;
		pik_hfs = NULL;
		partition_sfs = NULL;

		if (NEED_SFS) sfs = sfs_window(freqData, snps, SFS_SUB, CONST_N);

		int s = 0;
		int s_pi = MISSING; //note storage location of pi if it exists
		//useful for calculating Taj's D or F&W's H
		int s_S = MISSING;
		for (int j = 0; j < NOPTS; j++) {
			if (STATS[j].compare(ARG_PI) == 0 && params->getBoolFlag(ARG_PI)) {
				results[i][s] = pi_from_sfs(sfs);
				s_pi = s;
				s++;
			}
			else if (STATS[j].compare(ARG_PIK) == 0 && PIK_CHOICE[0] != 0) {
				pik_hfs = hfs_window(hapData, snps);
				for (int k = 0; k < PIK_CHOICE.size(); k++) {
					results[i][s] = pi_k2(pik_hfs, PIK_CHOICE[k]);
					s++;
				}
			}
			else if (STATS[j].compare(ARG_SEGSITES) == 0 && params->getBoolFlag(ARG_SEGSITES)) {
				results[i][s] = segsites(sfs);
				s_S = s;
				s++;
			}
			else if (STATS[j].compare(ARG_EHH) == 0 && EHH_WINS[0] != 0) {
				vector< pair_t* > *ehh_windows = getEHHWindows(snps->start, snps->winStart, WINSIZE, EHH_WINS, mapData, USE_BP);
				for (int w = 0; w < ehh_windows->size(); w++) {
					hfs = hfs_window(hapData, ehh_windows->at(w));
					results[i][s] = ehh_from_hfs(hfs);
					s++;
					releaseHaplotypeFrequencySpectrum(hfs);
				}
				releaseAllWindows(ehh_windows);
			}
			else if (STATS[j].compare(ARG_EHHK) == 0 && EHHK_CHOICES[0] != 0) {
				vector< pair_t* > *ehh_windows = getEHHWindows(snps->start, snps->winStart, WINSIZE, EHH_WINS, mapData, USE_BP);
				for (int w = 0; w < ehh_windows->size(); w++) {
					hfs = hfs_window(hapData, ehh_windows->at(w));
					for (int k = 0; k < EHHK_CHOICES.size(); k++) {
						results[i][s] = ehhk_from_hfs(hfs, EHHK_CHOICES[k]);
						s++;
					}
					releaseHaplotypeFrequencySpectrum(hfs);
				}
				releaseAllWindows(ehh_windows);
			}
			else if (STATS[j].compare(ARG_EHH_CM) == 0 && EHH_CM[0] != 0) {
				//subwindows measured in map units, centred on the parent
				//window's genetic midpoint
				vector< pair_t* > *cm_windows = getEHHWindowsGenetic(snps, EHH_CM, mapData);
				for (unsigned int w = 0; w < cm_windows->size(); w++) {
					hfs = hfs_window(hapData, cm_windows->at(w));
					results[i][s] = ehh_from_hfs(hfs);
					s++;
					releaseHaplotypeFrequencySpectrum(hfs);
				}
				releaseAllWindows(cm_windows);
			}
			else if (STATS[j].compare(ARG_TAJ_D) == 0 && params->getBoolFlag(ARG_TAJ_D)) {
				if (s_pi >= 0 && s_S >= 0) results[i][s] = tajimaD_from_sfs(sfs, results[i][s_pi], results[i][s_S]);
				else if (s_pi < 0 && s_S >= 0) results[i][s] = tajimaD_from_sfs(sfs, s_pi, results[i][s_S]);
				else if (s_pi >= 0 && s_S < 0) results[i][s] = tajimaD_from_sfs(sfs, results[i][s_pi], s_S);
				else results[i][s] = tajimaD_from_sfs(sfs);
				s++;
			}
			else if (STATS[j].compare(ARG_FAY_WU_H) == 0 && params->getBoolFlag(ARG_FAY_WU_H)) {
				if (s_pi >= 0) results[i][s] = fayWuH_from_sfs(sfs, results[i][s_pi]);
				else results[i][s] = fayWuH_from_sfs(sfs);
				s++;
			}
		}

		releaseArray(sfs);

		if (DO_PARTITION) {
			char part[2];
			part[0] = 'A';
			part[1] = '\0';
			vector< pair_t* > *partition_windows = getPartitionWindows(snps->start, snps->winStart, PARTITIONS, mapData, USE_BP);
			for (int p = 0; p < partition_windows->size(); p++) {
				int s_pi0 = MISSING;
				int s_S0 = MISSING;
				string partStr(part);
				partition_snps = partition_windows->at(p);
				partition_sfs = NULL;
				if (NEED_SFS) partition_sfs = sfs_window(freqData, partition_snps, SFS_SUB, CONST_N);

				for (int j = 0; j < NOPTS; j++) {
					if (STATS[j].compare(ARG_PI) == 0 && params->getBoolFlag(ARG_PI)) {
						results[i][s] = pi_from_sfs(partition_sfs);
						s_pi0 = s;
						s++;
					}
					else if (STATS[j].compare(ARG_PIK) == 0 && PIK_CHOICE[0] != 0) {
						pair_t *shifted_snps = new pair_t;
						shifted_snps->start = partition_snps->start - snps->start;
						shifted_snps->end = partition_snps->end - snps->start;
						for (int k = 0; k < PIK_CHOICE.size(); k++) {
							results[i][s] = pi_k2(pik_hfs, PIK_CHOICE[k], shifted_snps);
							s++;
						}
						delete shifted_snps;
					}
					else if (STATS[j].compare(ARG_SEGSITES) == 0 && params->getBoolFlag(ARG_SEGSITES)) {
						results[i][s] = segsites(partition_sfs);
						s_S0 = s;
						s++;
					}
					else if (STATS[j].compare(ARG_TAJ_D) == 0 && params->getBoolFlag(ARG_TAJ_D)) {
						if (s_pi0 >= 0 && s_S0 >= 0) results[i][s] = tajimaD_from_sfs(partition_sfs, results[i][s_pi0], results[i][s_S0]);
						else if (s_pi0 < 0 && s_S0 >= 0) results[i][s] = tajimaD_from_sfs(partition_sfs, s_pi0, results[i][s_S0]);
						else if (s_pi0 >= 0 && s_S0 < 0) results[i][s] = tajimaD_from_sfs(partition_sfs, results[i][s_pi0], s_S0);
						else results[i][s] = tajimaD_from_sfs(partition_sfs);
						s++;
					}
					else if (STATS[j].compare(ARG_FAY_WU_H) == 0 && params->getBoolFlag(ARG_FAY_WU_H)) {
						if (s_pi0 >= 0) results[i][s] = fayWuH_from_sfs(partition_sfs, results[i][s_pi0]);
						else results[i][s] = fayWuH_from_sfs(partition_sfs);
						s++;
					}
					else if (STATS[j].compare(ARG_EHH) == 0 && EHH_WINS[0] != 0 && params->getBoolFlag(ARG_EHH_PART)) {
						hfs = hfs_window(hapData, partition_snps);
						results[i][s] = ehh_from_hfs(hfs);
						s++;
						releaseHaplotypeFrequencySpectrum(hfs);
					}
					else if (STATS[j].compare(ARG_EHHK) == 0 && EHHK_CHOICES[0] != 0 && params->getBoolFlag(ARG_EHH_PART)) {

						hfs = hfs_window(hapData, partition_snps);
						for (int k = 0; k < EHHK_CHOICES.size(); k++) {
							results[i][s] = ehhk_from_hfs(hfs, EHHK_CHOICES[k]);
							s++;
						}
						releaseHaplotypeFrequencySpectrum(hfs);

					}
				}
				part[0]++;
				releaseArray(partition_sfs);
			}
		}
		releaseHaplotypeFrequencySpectrum(pik_hfs);

		//The header is built independently in main, by buildColumnNames, from
		//the same flags. If the two ever disagree the columns silently
		//misalign with their names, so check rather than assume.
		if (s != numStats) {
			cerr << "ERROR: window " << i << " produced " << s
			     << " statistics but the header declares " << numStats
			     << "; buildColumnNames and calc_stats have drifted apart.\n";
			throw 0;
		}
	}
	return;
}

//EHH subwindows placed by GENETIC distance. The parent window's genetic
//midpoint is the mean of its first and last SNP's map positions, and each
//subwindow takes the SNPs within +/- width/2 of that point.
//
//This is the only thing in divstats that reads geneticPos. Before --ehh-cm,
//a genetic map was mandatory for every EHH calculation and then never
//consulted: subwindows were placed by physical position in both --sites and
//--bp mode, so users produced a recombination map for a computation that
//ignored it. Widths in map units are comparable between regions of differing
//recombination rate, which a width in base pairs is not.
vector< pair_t* > *getEHHWindowsGenetic(pair_t *parentWin, vector<double> &EHH_CM, MapData *mapData) {
	vector< pair_t* > *ehh_windows = new vector< pair_t* >;
	int numSnps = mapData->nloci;

	for (unsigned int i = 0; i < EHH_CM.size(); i++) {
		pair_t *snps = new pair_t;

		//An empty parent window has no midpoint to centre on.
		if (parentWin->end < parentWin->start) {
			snps->start = parentWin->start;
			snps->end = parentWin->start - 1;
			snps->winStart = parentWin->winStart;
			snps->winEnd = parentWin->winEnd;
			ehh_windows->push_back(snps);
			continue;
		}

		double mid = 0.5 * (mapData->geneticPos[parentWin->start] +
		                    mapData->geneticPos[parentWin->end]);
		double lo = mid - 0.5 * EHH_CM[i];
		double hi = mid + 0.5 * EHH_CM[i];

		//Clamp to the parent window: a subwindow may not reach outside it,
		//matching the constraint --ehh already documents.
		int a = parentWin->start;
		while (a <= parentWin->end && mapData->geneticPos[a] < lo) a++;
		int b = a;
		while (b <= parentWin->end && b < numSnps && mapData->geneticPos[b] <= hi) b++;
		b--;

		snps->start = a;
		snps->end = b;	//end < start means no SNP falls in the subwindow
		snps->winStart = (b >= a) ? mapData->physicalPos[a] : parentWin->winStart;
		snps->winEnd   = (b >= a) ? mapData->physicalPos[b] : parentWin->winEnd;
		ehh_windows->push_back(snps);
	}
	return ehh_windows;
}

//%g gives the shortest round-trip-ish form, so 0.05 stays "0.05" and 1.0
//becomes "1" -- stable column names without trailing zeros.
string dbl2str(double d) {
	char buf[32];
	snprintf(buf, sizeof(buf), "%g", d);
	return string(buf);
}

string int2str(int i) {
	char buffer[10];
	sprintf(buffer, "%d", i);
	return string(buffer);
}

vector< pair_t* > *getPartitionWindows(int snpStart, int winStart, vector<int> &PARTITIONS, MapData *mapData, bool USE_BP) {
	vector< pair_t* > *partition_windows = new vector< pair_t* >;
	if (USE_BP){
		int partitionSnpIndexStart = snpStart;
		int partitionCurrWinStart = winStart;
		for (int i = 0; i < PARTITIONS.size(); i++) {
			pair_t *partition_snps = findInclusiveSNPIndicies(partitionSnpIndexStart, partitionCurrWinStart, PARTITIONS[i], mapData);
			partition_windows->push_back(partition_snps);
			partitionSnpIndexStart = partition_snps->end;
			partitionCurrWinStart += PARTITIONS[i];
		}
	}
	else{//USE_SITES
		int numSnps = mapData->nloci;
		int currStart = snpStart;
		int currEnd = -1;
		for (int i = 0; i < PARTITIONS.size(); i++){
			pair_t *partition_snps = new pair_t;
			currEnd = currStart + PARTITIONS[i] - 1;
			partition_snps->start = currStart;
			partition_snps->end = (currEnd >= numSnps) ? numSnps -1 : currEnd;
			partition_windows->push_back(partition_snps);
			currStart = currEnd + 1;
		}
	}
	return partition_windows;
}

vector< pair_t* > *getEHHWindows(int snpStart, int winStart, int WINSIZE, vector<int> &EHH_WINS, MapData *mapData, bool USE_BP) {
	vector< pair_t* > *ehh_windows = new vector< pair_t* >;
	if(USE_BP){
		int currWinStart = winStart;
		for (int i = 0; i < EHH_WINS.size(); i++) {
			pair_t *snps = findInclusiveSNPIndicies(snpStart, ( winStart + (WINSIZE * 0.5) - (EHH_WINS[i] * 0.5) ) , EHH_WINS[i], mapData);
			ehh_windows->push_back(snps);
		}
	}
	else{//USE_SITES
		double mid = (WINSIZE - 1) * 0.5 + snpStart;
		for (int i = 0; i < EHH_WINS.size(); i++) {
			pair_t *snps = new pair_t;
			if (EHH_WINS[i] % 2 == 0){
				snps->start = int(mid - (EHH_WINS[i] * 0.5) + 0.5);
				snps->end = int(mid + (EHH_WINS[i] * 0.5));
			}
			else{
				if (WINSIZE % 2 == 0){
					snps->start = int(mid - (EHH_WINS[i] * 0.5) + 1);
					snps->end = int(mid + (EHH_WINS[i] * 0.5) - 1);
				}
				else{
					snps->start = int(mid - (EHH_WINS[i] * 0.5) + 0.5);
					snps->end = int(mid + (EHH_WINS[i] * 0.5));
				}
			}
			ehh_windows->push_back(snps);
		}
	}
	return ehh_windows;
}


pair_t* findInclusiveSNPIndicies(int startSnpIndex, int currWinStart, int WINSIZE, MapData* mapData) {

	int currWinEnd = currWinStart + WINSIZE - 1;
	int endSnpIndex = startSnpIndex;
	int numSnps = mapData->nloci;

	pair_t* snps = new pair_t;
	snps->winStart = currWinStart;
	snps->winEnd = currWinEnd;
	if (mapData->physicalPos[numSnps - 1] < currWinStart) {
		snps->start = numSnps;
		snps->end = numSnps - 1;
		return snps;
	}

	//The window is [currWinStart, currWinEnd] inclusive at both ends, as the
	//function name says. The previous walk was
	//
	//	while (physicalPos[endSnpIndex] < currWinEnd) { endSnpIndex++; ... }
	//	endSnpIndex--;
	//
	//which dropped a SNP sitting exactly on currWinEnd: the loop stops as soon
	//as the position is not strictly less than the end, so for equality it
	//stops ON the SNP that belongs in the window, and the decrement then
	//excludes it. With --bp and a round --winsize, positions landing exactly on
	//a boundary are common, and such a SNP was counted in neither the window
	//that should have contained it nor the next one.
	//
	//The bound test also broke the last window. It read
	//`if (endSnpIndex >= numSnps - 1) break;` after the increment, so the walk
	//stopped at numSnps - 1 and the unconditional decrement then gave
	//numSnps - 2 -- discarding the final SNP of the data whenever the last
	//window extended past it. The clamp on the next line tested
	//`endSnpIndex >= numSnps`, which that break made unreachable.
	while (startSnpIndex < numSnps && mapData->physicalPos[startSnpIndex] < currWinStart) {
		startSnpIndex++;
	}
	endSnpIndex = startSnpIndex;
	while (endSnpIndex < numSnps && mapData->physicalPos[endSnpIndex] <= currWinEnd) {
		endSnpIndex++;
	}
	endSnpIndex--;	//last SNP at or before currWinEnd; startSnpIndex-1 if none

	snps->start = startSnpIndex;
	snps->end = endSnpIndex;	//end < start means the window holds no SNPs
	return snps;
}