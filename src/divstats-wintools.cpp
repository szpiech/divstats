#include "divstats-wintools.h"

vector< pair_t* > *findAllWindows(MapData *mapData, int WINSIZE, int WINSTEP) {
	int currWinStart = 1;//mapData->physicalPos[0];
	int currWinEnd = currWinStart + WINSIZE - 1;
	int numSnps = mapData->nloci;
	int endOfData = mapData->physicalPos[numSnps - 1];
	pair_t *snpIndex = new pair_t;
	snpIndex->start = 0;
	snpIndex->end = -1;
	int numInWindow;

	vector< pair_t* > *windows = new vector< pair_t* >;
	for (currWinStart; currWinStart < endOfData; currWinStart += WINSTEP, currWinEnd += WINSTEP) {

		//Find SNP index boundaries for the whole window
		pair_t *snps = findInclusiveSNPIndicies(snpIndex->start, currWinStart, WINSIZE, mapData);
		windows->push_back(snps);

		snpIndex->start = snps->start;
		snpIndex->end = snps->end;
	}

	return windows;
}

void releaseAllWindows(vector< pair_t* > *windows) {
	for (int i = 0; i < windows->size(); i++) delete windows->at(i);
	delete windows;
	return;
}

void calc_stats(void *order) {
	work_order_t *p = (work_order_t *)order;
	HaplotypeData *hapData = p->hapData;
	MapData *mapData = p->mapData;
	FreqData *freqData = p->freqData;
	param_t *params = p->params;
	vector< pair_t* > *windows = p->windows;
	double **results = p->results;
	int id = p->id;
	vector<int> PIK_CHOICE = params->getIntListFlag(ARG_PIK);
	vector<int> EHH_WINS = params->getIntListFlag(ARG_EHH);
	vector<int> EHHK_CHOICES = params->getIntListFlag(ARG_EHHK);
	vector<int> PARTITIONS = params->getIntListFlag(ARG_PARTITION);
	bool DO_PARTITION = p->DO_PARTITION;

	int numThreads = params->getIntFlag(ARG_THREADS);
	array_t *sfs, *partition_sfs;
	HaplotypeFrequencySpectrum *hfs, *partition_hfs;
	pair_t *snps, *partition_snps;


	//Cycle over all windows and calculate stats
	for (int i = id; i < windows->size(); i += numThreads) {
		snps = windows->at(i);
		sfs = sfs_window(freqData, snps);

		int s = 0;
		for (int j = 0; j < NOPTS; j++) {
			if (STATS[j].compare(ARG_PI) == 0 && params->getBoolFlag(ARG_PI)) {
				results[i][s] = pi_from_sfs(sfs);
				s++;
			}
			else if (STATS[j].compare(ARG_PIK) == 0 && PIK_CHOICE[0] != 0) {
				hfs = hfs_window(hapData, snps);
				for (int k = 0; k < PIK_CHOICE.size(); k++) {
					results[i][s] = pi_k2(hfs, PIK_CHOICE[k]);
					s++;
				}
				releaseHaplotypeFrequencySpectrum(hfs);
			}
			else if (STATS[j].compare(ARG_SEGSITES) == 0 && params->getBoolFlag(ARG_SEGSITES)) {
				results[i][s] = segsites(sfs);
				s++;
			}
			else if (STATS[j].compare(ARG_EHH) == 0 && EHH_WINS[0] != 0) {
				for (int k = 0; k < EHH_WINS.size(); k++) {
					results[i][s] = -9;//ehh();
					s++;
				}
			}
			else if (STATS[j].compare(ARG_EHHK) == 0 && EHHK_CHOICES[0] != 0) {
				for (int k = 0; k < EHHK_CHOICES.size(); k++) {
					results[i][s] = -9;//ehhk();
					s++;
				}
			}
			else if (STATS[j].compare(ARG_TAJ_D) == 0 && params->getBoolFlag(ARG_TAJ_D)) {
				results[i][s] = -9;
				s++;
			}
			else if (STATS[j].compare(ARG_FAY_WU_H) == 0 && params->getBoolFlag(ARG_FAY_WU_H)) {
				results[i][s] = -9;
				s++;
			}

			if (DO_PARTITION) {
				vector< pair_t* > *partition_windows = getPartitionWindows(snps->start, snps->winStart, PARTITIONS, mapData);
				for (int p = 0; p < partition_windows->size(); p++) {
					partition_snps = partition_windows->at(p);
					partition_sfs = sfs_window(freqData, partition_snps);

					if (STATS[j].compare(ARG_PI) == 0 && params->getBoolFlag(ARG_PI)) {
						results[i][s] = pi_from_sfs(partition_sfs);
						s++;
					}
					else if (STATS[j].compare(ARG_PIK) == 0 && PIK_CHOICE[0] != 0) {
						partition_hfs = hfs_window(hapData, partition_snps);
						for (int k = 0; k < PIK_CHOICE.size(); k++) {
							results[i][s] = pi_k2(partition_hfs, PIK_CHOICE[k]);
							s++;
						}
						releaseHaplotypeFrequencySpectrum(partition_hfs);
					}
					else if (STATS[j].compare(ARG_SEGSITES) == 0 && params->getBoolFlag(ARG_SEGSITES)) {
						results[i][s] = segsites(partition_sfs);
						s++;
					}
					else if (STATS[j].compare(ARG_EHH) == 0 && EHH_WINS[0] != 0) {
						for (int k = 0; k < EHH_WINS.size(); k++) {
							results[i][s] = -9;//ehh();
							s++;
						}
					}
					else if (STATS[j].compare(ARG_EHHK) == 0 && EHHK_CHOICES[0] != 0) {
						for (int k = 0; k < EHHK_CHOICES.size(); k++) {
							results[i][s] = -9;//ehhk();
							s++;
						}
					}
					else if (STATS[j].compare(ARG_TAJ_D) == 0 && params->getBoolFlag(ARG_TAJ_D)) {
						results[i][s] = -9;
						s++;
					}
					else if (STATS[j].compare(ARG_FAY_WU_H) == 0 && params->getBoolFlag(ARG_FAY_WU_H)) {
						results[i][s] = -9;
						s++;
					}

					releaseArray(partition_sfs);
				}
			}
		}
		releaseArray(sfs);
	}
	return;
}

vector< pair_t* > *getPartitionWindows(int snpStart, int winStart, vector<int> &PARTITIONS, MapData *mapData) {
	vector< pair_t* > *partition_windows = new vector< pair_t* >;
	pair_t *partitionSnpIndex = new pair_t;
	partitionSnpIndex->start = snpStart;
	partitionSnpIndex->end = snpStart - 1;
	int partitionCurrWinStart = winStart;
	for (int i = 0; i < PARTITIONS.size(); i++) {
		pair_t *partition_snps = findInclusiveSNPIndicies(partitionSnpIndex->start, partitionCurrWinStart, PARTITIONS[i], mapData);
		partition_windows->push_back(partition_snps);
		partitionSnpIndex->start = partition_snps->start;
		partitionSnpIndex->end = partition_snps->end;
		partitionCurrWinStart += PARTITIONS[i];
	}
	delete partitionSnpIndex;
	return partition_windows;
}

pair_t* findInclusiveSNPIndicies(int startSnpIndex, int currWinStart, int WINSIZE, MapData* mapData) {

	int currWinEnd = currWinStart + WINSIZE - 1;
	int endSnpIndex = startSnpIndex;
	int numSnps = mapData->nloci;

	pair_t* snps = new pair_t;
	snps->winStart = currWinStart;
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