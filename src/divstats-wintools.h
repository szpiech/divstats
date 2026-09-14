#ifndef __DIVSTATS_WINTOOLS_H__
#define __DIVSTATS_WINTOOLS_H__

#include <vector>
#include <pthread.h>
#include <ctime>
#include <map>
#include "divstats-data.h"
#include "divstats-winstats.h"
#include "param_t.h"
#include "divstats-cli.h"

//Shared progress counter for the worker threads.
struct progress_t
{
    pthread_mutex_t lock;
    long done;
    long total;
    int lastTenth;       //highest 10% mark already reported
    time_t start;
    bool announced;      //true once anything has been printed
};

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

    //Progress reporting. The window loop used to print nothing at all
    //between "Calculating N statistics in M windows." and the output file,
    //so a multi-hour run gave no sign of life. Threads share one counter
    //under a mutex -- the lock is taken once per window, which is nothing
    //next to the window's own work, and it keeps the count exact rather
    //than reading another thread's int without synchronisation.
    progress_t *progress;

    param_t *params;

    bool DO_PARTITION;
    bool USE_BP;
    bool SFS_SUB;
    int TARGET_N;            //0 = each window's own minimum

    //per-window metadata, filled by calc_stats and written by main
    int *nhapsUsed;
    int *nSitesUsed;
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
string dbl2str(double d);

//EHH subwindows placed by GENETIC distance: width is in map units, centred on
//the genetic midpoint of the parent window. Requires a --map.
vector< pair_t* > *getEHHWindowsGenetic(pair_t *parentWin, vector<double> &EHH_CM, MapData *mapData);

#endif