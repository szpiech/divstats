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
    string *names;
    //ofstream *flog;
    //Bar *bar;

    param_t *params;

    bool DO_PARTITION;
};

pair_t* findInclusiveSNPIndicies(int startSnpIndex, int currWinStart, int WINSIZE, MapData* mapData);

vector< pair_t* > *findAllWindows(MapData *mapData, int WINSIZE, int WINSTEP);
void releaseAllWindows(vector< pair_t* > *windows);

vector< pair_t* > *getPartitionWindows(int snpStart, int winStart, vector<int> &PARTITIONS, MapData *mapData);

vector< pair_t* > *getEHHWindows(int snpStart, int winStart, int WINSIZE, vector<int> &EHH_WINS, MapData *mapData);

void calc_stats(void *work_order);

string int2str(int i);

#endif