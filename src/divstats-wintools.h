#ifndef __DIVSTATS_WINTOOLS_H__
#define __DIVSTATS_WINTOOLS_H__

#include <vector>
#include <map>
#include "divstats-data.h"
#include "divstats-winstats.h"
#include "param_t.h"
#include "divstats-cli.h"

struct work_order_t
{
    int id;
    int numStats;
    HaplotypeData *hapData;
    MapData *mapData;
    FreqData *freqData;

	vector< pair_t* > *windows;

    double **results;
    //ofstream *flog;
    //Bar *bar;

    param_t *params;

    bool DO_PARTITION;
    bool USE_BP;
    bool SFS_SUB;
};

pair_t* findInclusiveSNPIndicies(int startSnpIndex, int currWinStart, int WINSIZE, MapData* mapData);

vector< pair_t* > *findAllWindows(MapData *mapData, int WINSIZE, int WINSTEP, bool USE_BP);
void releaseAllWindows(vector< pair_t* > *windows);

vector< pair_t* > *getPartitionWindows(int snpStart, int winStart, vector<int> &PARTITIONS, MapData *mapData, bool USE_BP);

vector< pair_t* > *getEHHWindows(int snpStart, int winStart, int WINSIZE, vector<int> &EHH_WINS, MapData *mapData, bool USE_BP);

//Output column names, in the order calc_stats fills results[][]. A pure
//function of the command line, so main can write the header before any
//thread starts.
vector<string> buildColumnNames(param_t *params, bool DO_PARTITION);

void calc_stats(void *work_order);

string int2str(int i);

#endif