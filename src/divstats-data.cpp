/* selscan -- a program to calculate EHH-based scans for positive selection in genomes
   Copyright (C) 2014  Zachary A Szpiech

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
#include "divstats-data.h"

HaplotypeFrequencySpectrum *initHaplotypeFrequencySpectrum(){
    HaplotypeFrequencySpectrum *hfs = new HaplotypeFrequencySpectrum;
    hfs->sortedCount = NULL;
    hfs->size = 0;
    hfs->numUniq = 0;
    return hfs;
}

void releaseHaplotypeFrequencySpectrum(HaplotypeFrequencySpectrum *hfs){
    if(hfs == NULL){
        return;
    }

    if(hfs->sortedCount != NULL){
        delete [] hfs->sortedCount;
    }

    delete hfs;
    return;
}

array_t *initArray(int size, double fill){
    array_t *data = new array_t;
    data->size = size;
    data->data = new double[size];
    for(int i = 0; i < size; i++){
        data->data[i] = fill;
    }
    return data;
}
void releaseArray(array_t* data){
    if(data == NULL){
        return;
    }

    if(data->data != NULL){
        delete [] data->data;
    }
    delete data;
    return;
}

FreqData *initFreqData(int nhaps, int nloci) {
    FreqData *freqData = new FreqData;
    freqData->count = new int[nloci];
    freqData->nmissing = new int[nloci];
    for (int i = 0; i < nloci; i++) {
        freqData->count[i] = MISSING;
    }
    freqData->nloci = nloci;
    freqData->nhaps = nhaps;
    freqData->maxMissing = -1;
    return freqData;
}

FreqData *initFreqData(HaplotypeData* data) {
    if (data == NULL) {
        cerr << "ERROR: Can not compute frequencies on NULL data.\n";
        throw 0;
    }
    FreqData *freqData = initFreqData(data->nhaps,data->nloci);
    //cerr << "Calculating frequencies on " << data->nhaps << " haps across " << data->nloci << " loci.\n";

    for (int locus = 0; locus < data->nloci; locus++)
    {
        freqData->count[locus] = 0;
        freqData->nmissing[locus] = 0;
        for (int hap = 0; hap < data->nhaps; hap++)
        {
            char a = hapGet(data->data[hap], locus);
            freqData->count[locus] += ( a == '1' ? 1 : 0 );
            freqData->nmissing[locus] += ( a == MISSING_ALLELE ? 1 : 0 );
            if (a != '0' && a != '1' && a != MISSING_ALLELE)
            {
                cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                throw 0;
            }
        }
        if (freqData->maxMissing < freqData->nmissing[locus]) freqData->maxMissing = freqData->nmissing[locus];
    }

    return freqData;
}

void releaseFreqData(FreqData *data){
    if(data == NULL){
        return;
    }
    if(data->count != NULL){
        delete [] data->count;
    }
    if(data->nmissing != NULL){
        delete [] data->nmissing;
    }
    delete data;
    return;
}


//reads in map data and also does basic checks on integrity of format
//returns a populated MapData structure if successful
//throws an exception otherwise
MapData *readMapData(string filename, int expected_loci)
{
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    //int fileStart = fin.tellg();
    string line;
    int nloci = 0;
    int num_cols = 4;
    int current_cols = 0;
    while (getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);
        if (current_cols != num_cols)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_cols
                 << ", but expected " << num_cols << ".\n";
            throw 0;
        }
    }

    if (nloci != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci << ".\n";
        throw 0;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading map data for " << nloci << " loci\n";

    MapData *data = initMapData(nloci);

    string chr;
    for (int locus = 0; locus < data->nloci; locus++)
    {
        fin >> data->chr;
        fin >> data->locusName[locus];
        fin >> data->geneticPos[locus];
        fin >> data->physicalPos[locus];
        //if (USE_PMAP) data->geneticPos[locus] = data->physicalPos[locus];
    }

    fin.close();
    return data;
}

MapData *readMapDataTPED(string filename, int expected_loci, int expected_haps)
{
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    //int fileStart = fin.tellg();
    string line;
    int nloci = 0;
    int num_cols = 4;
    int current_cols = 0;
    while (getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);
        if (current_cols != num_cols + expected_haps)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_cols
                 << ", but expected " << num_cols + expected_haps << ".\n";
            throw 0;
        }
    }

    if (nloci != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci << ".\n";
        throw 0;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading map data for " << nloci << " loci\n";

    MapData *data = initMapData(nloci);

    string chr;
    for (int locus = 0; locus < data->nloci; locus++)
    {
        fin >> data->chr;
        fin >> data->locusName[locus];
        fin >> data->geneticPos[locus];
        fin >> data->physicalPos[locus];
        getline(fin, line);
    }

    fin.close();
    return data;
}

//allocates the arrays and populates them with MISSING or "--" depending on type
MapData *initMapData(int nloci)
{
    if (nloci < 1)
    {
        cerr << "ERROR: number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }

    MapData *data = new MapData;
    data->nloci = nloci;
    data->locusName = new string[nloci];
    data->physicalPos = new int[nloci];
    data->geneticPos = new double[nloci];

    for (int locus = 0; locus < nloci; locus++)
    {
        data->locusName[locus] = "--";
        data->physicalPos[locus] = MISSING;
        data->geneticPos[locus] = MISSING;
    }

    return data;
}

void releaseMapData(MapData *data)
{
    if (data == NULL) return;
    data->nloci = -9;
    delete [] data->locusName;
    delete [] data->physicalPos;
    delete [] data->geneticPos;
    delete data;
    data = NULL;
    return;
}

HaplotypeData *readHaplotypeDataTPED(string filename)
{
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int numMapCols = 4;
    //int fileStart = fin.tellg();
    string line;
    int nloci = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;
    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        //getline(fin,line);
        //if(fin.eof()) break;
        nloci++;
        current_nhaps = countFields(line);
        //cout << "nloci: " << current_nhaps << endl;
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_nhaps
                 << " fields, but the previous line has " << previous_nhaps << " fields.\n";
            throw 0;
        }
        previous_nhaps = current_nhaps;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading " << current_nhaps - numMapCols << " haplotypes and " << nloci << " loci...\n";

    HaplotypeData *data = initHaplotypeData(current_nhaps - numMapCols, nloci);

    string junk;
    string allele;
    for (int locus = 0; locus < data->nloci; locus++)
    {
        for (int i = 0; i < numMapCols; i++)
        {
            fin >> junk;
        }
        for (int hap = 0; hap < data->nhaps; hap++)
        {
            fin >> allele;
            if (allele.compare("0") != 0 && allele.compare("1") != 0 && allele.compare(TPED_MISSING) != 0)
            {
                cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                throw 0;
            }
            if (allele.compare(TPED_MISSING) == 0) hapSet(data->data[hap], locus, MISSING_ALLELE);
            else hapSet(data->data[hap], locus, allele[0]);
        }
    }

    fin.close();

    return data;
}


HaplotypeData *initHaplotypeData(unsigned int nhaps, unsigned int nloci)
{
    if (nhaps < 1 || nloci < 1)
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }

    HaplotypeData *data = new HaplotypeData;
    data->nhaps = nhaps;
    data->nloci = nloci;

    data->data = new char *[nhaps];
    for (unsigned int i = 0; i < nhaps; i++)
    {
        //packed 2 bits per site; the fill is MISSING_CHAR in all four slots
        long nbytes = hapRowBytes((long)nloci);
        data->data[i] = new char[nbytes];
        for (long j = 0; j < nbytes; j++)
        {
            data->data[i][j] = (char)0xFF;      //four copies of code 3 = '9'
        }
    }

    return data;
}

void releaseHapData(HaplotypeData *data)
{
    if (data == NULL) return;
    for (int i = 0; i < data->nhaps; i++)
    {
        delete [] data->data[i];
    }

    delete [] data->data;

    //The struct itself was never freed. The three assignments below it -- and
    //`data = NULL`, which only clears this function's own by-value copy of the
    //caller's pointer and so does nothing at all -- left the HaplotypeData
    //allocation leaked on every call. releaseMapData, releaseFreqData and
    //releaseArray all delete their struct; this one did not.
    delete data;
    return;
}


int countFields(const string &str)
{
    string::const_iterator it;
    int result;
    int numFields = 0;
    int seenChar = 0;
    for (it = str.begin() ; it < str.end(); it++)
    {
        result = isspace(*it);
        if (result == 0 && seenChar == 0)
        {
            numFields++;
            seenChar = 1;
        }
        else if (result != 0)
        {
            seenChar = 0;
        }
    }
    return numFields;
}



