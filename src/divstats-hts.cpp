/* divstats -- htslib-backed variant reader
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
#include "divstats-hts.h"

#include <vector>
#include <sstream>

#include "htslib/vcf.h"
#include "htslib/hts.h"

using namespace std;

//Small RAII holder so every error path closes the file and frees htslib's
//allocations. The reader throws 0 in a dozen places; without this each one
//would need its own cleanup.
namespace {
struct HtsReader {
   htsFile *fp;
   bcf_hdr_t *hdr;
   bcf1_t *rec;
   int32_t *gt;
   int ngt_arr;
   HtsReader() : fp(NULL), hdr(NULL), rec(NULL), gt(NULL), ngt_arr(0) {}
   ~HtsReader() {
      if (gt) free(gt);
      if (rec) bcf_destroy(rec);
      if (hdr) bcf_hdr_destroy(hdr);
      if (fp) hts_close(fp);
   }
};
}

void readVariantDataHTS(string filename, bool HEMI, int nThreads,
                        HaplotypeData **hapDataOut, MapData **mapDataOut)
{
   HtsReader R;

   cerr << "Opening " << filename << "...\n";
   R.fp = hts_open(filename.c_str(), "r");
   if (R.fp == NULL) {
      cerr << "ERROR: Failed to open " << filename << " for reading.\n";
      throw 0;
   }

   const htsFormat *fmt = hts_get_format(R.fp);
   if (fmt->category != variant_data) {
      cerr << "ERROR: " << filename << " is not variant data ("
           << hts_format_description(fmt) << ").\n";
      throw 0;
   }

   //Decoding runs in htslib's own thread pool when asked. This parallelises
   //BGZF block inflation, so it does something for a bgzipped VCF or a BCF
   //and nothing for plain gzip, which is a single deflate stream.
   if (nThreads > 1) hts_set_opt(R.fp, HTS_OPT_NTHREADS, nThreads);

   R.hdr = bcf_hdr_read(R.fp);
   if (R.hdr == NULL) {
      cerr << "ERROR: Failed to read a VCF/BCF header from " << filename << ".\n";
      throw 0;
   }

   int nsmpl = bcf_hdr_nsamples(R.hdr);
   if (nsmpl <= 0) {
      cerr << "ERROR: " << filename << " contains no samples.\n";
      throw 0;
   }

   int nhaps = HEMI ? nsmpl : 2 * nsmpl;

   //Genotypes accumulate as the file is read, because the locus count is not
   //known until the end. The old reader learned it by reading the whole file
   //first; growing this is what removes that pass.
   //
   //Accumulated SITE-MAJOR, in one buffer, and transposed once at the end.
   //Appending per haplotype instead -- vector<vector<char> > haps(nhaps), one
   //push_back per haplotype per record -- touched nhaps separate buffers at
   //every record: 20 million scattered appends on a 400-haplotype file, each
   //landing on a different cache line, plus nhaps independent growth
   //schedules. Site-major makes each record one contiguous run of nhaps
   //bytes, and the transpose is done in cache-sized blocks.
   vector<char> siteMajor;
   siteMajor.reserve((size_t)nhaps * 8192);
   vector<int> pos;
   vector<string> ids;

   R.rec = bcf_init();
   int rid = -1;              //chromosome of the first record
   string chrName;
   long nrec = 0;
   int lastPos = -1;

   while (bcf_read(R.fp, R.hdr, R.rec) == 0) {
      bcf_unpack(R.rec, BCF_UN_STR);
      nrec++;

      //--- B8: one chromosome per run -------------------------------------
      //MapData carries a single chr string, and the old reader simply
      //overwrote it on every record. Records from different chromosomes were
      //concatenated into one coordinate space, windows straddled the boundary,
      //and every output row was labelled with the last chromosome seen -- with
      //no warning. Refuse the input instead.
      if (rid < 0) {
         rid = R.rec->rid;
         chrName = bcf_hdr_id2name(R.hdr, R.rec->rid);
      }
      else if (R.rec->rid != rid) {
         cerr << "ERROR: " << filename << " contains more than one chromosome ("
              << chrName << " then " << bcf_hdr_id2name(R.hdr, R.rec->rid)
              << ", at record " << nrec << ").\n"
              << "       divstats analyses one chromosome per run: windows are\n"
              << "       placed in a single coordinate space and the output\n"
              << "       carries one chromosome label. Split the input first,\n"
              << "       e.g. bcftools view -r " << chrName << ".\n";
         throw 0;
      }

      int p = (int)(R.rec->pos + 1);       //htslib is 0-based, VCF is 1-based

      //--- B9: positions must be ascending --------------------------------
      //findInclusiveSNPIndicies walks physicalPos monotonically and assumes it
      //is sorted; nothing checked. An unsorted file was accepted silently, the
      //first window reported nSNPs = 0 with -999 statistics, and only a subset
      //of sites was ever visited.
      if (p <= lastPos) {
         cerr << "ERROR: " << filename << " is not sorted by position: record "
              << nrec << " is at " << chrName << ":" << p
              << ", after " << chrName << ":" << lastPos << ".\n"
              << "       Sort it first, e.g. bcftools sort.\n";
         throw 0;
      }
      lastPos = p;

      if (R.rec->n_allele > 2) {
         cerr << "ERROR: " << chrName << ":" << p << " has " << R.rec->n_allele
              << " alleles. divstats requires biallelic sites; filter first,\n"
              << "       e.g. bcftools view -m2 -M2.\n";
         throw 0;
      }

      //bcf_get_genotypes finds GT wherever it sits in FORMAT. The old reader
      //took characters 0 and 2 of each sample field without consulting FORMAT
      //at all, so a file with FORMAT=DP:GT failed with the misleading message
      //"Alleles must be coded 0/1 or missing only".
      int ngt = bcf_get_genotypes(R.hdr, R.rec, &R.gt, &R.ngt_arr);
      if (ngt <= 0) {
         cerr << "ERROR: no GT field at " << chrName << ":" << p
              << " (FORMAT has no GT).\n";
         throw 0;
      }
      int maxPloidy = ngt / nsmpl;
      int wanted = HEMI ? 1 : 2;
      if (maxPloidy < wanted) {
         cerr << "ERROR: " << chrName << ":" << p << " has ploidy " << maxPloidy
              << " but " << (HEMI ? "1 haplotype" : "2 haplotypes")
              << " per sample is required"
              << (HEMI ? ".\n" : "; use --hemi for haploid data.\n");
         throw 0;
      }

      //One resize per record, then plain stores: push_back per allele paid a
      //capacity check on each of the nhaps writes.
      size_t base = siteMajor.size();
      siteMajor.resize(base + (size_t)nhaps);
      char *w = &siteMajor[base];
      int wi = 0;

      for (int i = 0; i < nsmpl; i++) {
         int32_t *g = R.gt + i * maxPloidy;
         for (int j = 0; j < wanted; j++) {
            char allele;
            if (g[j] == bcf_int32_vector_end || bcf_gt_is_missing(g[j])) {
               allele = MISSING_ALLELE;
            }
            else {
               int idx = bcf_gt_allele(g[j]);
               if (idx < 0 || idx > 1) {
                  cerr << "ERROR: allele index " << idx << " at " << chrName
                       << ":" << p << "; alleles must be coded 0/1 or missing.\n";
                  throw 0;
               }
               allele = (char)('0' + idx);
            }
            w[wi++] = allele;
         }
      }

      pos.push_back(p);
      const char *id = R.rec->d.id;
      ids.push_back((id && id[0]) ? string(id) : string("."));
   }

   long nloci = (long)pos.size();
   if (nloci <= 0) {
      cerr << "ERROR: " << filename << " contains no variant records.\n";
      throw 0;
   }

   cerr << "Loading " << nhaps << " haplotypes and " << nloci << " loci...\n";

   HaplotypeData *hapData = initHaplotypeData(nhaps, (unsigned int)nloci);

   //Transpose site-major -> haplotype-major. Done in blocks of loci so that
   //the source slab being read (BLOCK * nhaps bytes) stays resident while all
   //nhaps destination rows are filled from it; the destination writes are
   //contiguous within a row. A straight nhaps x nloci transpose instead
   //re-reads the whole source once per haplotype.
   //Blocks are a multiple of 4 loci so that a block boundary never falls
   //inside a packed byte; each destination byte is then written once, by one
   //block, and the four sites feeding it are all in this block.
   const long BLOCK = 256;
   for (long b = 0; b < nloci; b += BLOCK) {
      long bend = (b + BLOCK < nloci) ? b + BLOCK : nloci;
      for (int h = 0; h < nhaps; h++) {
         char *dst = hapData->data[h];
         const char *src = &siteMajor[(size_t)b * nhaps + h];
         for (long l = b; l < bend; l++) {
            hapSet(dst, l, *src);
            src += nhaps;
         }
      }
   }

   MapData *mapData = initMapData((int)nloci);
   mapData->chr = chrName;
   for (long i = 0; i < nloci; i++) {
      mapData->physicalPos[i] = pos[i];
      mapData->locusName[i]   = ids[i];
      //matches readMapDataVCF: with no genetic map, cM is set to bp. Nothing
      //in the tree reads geneticPos, but keep the behaviour identical.
      mapData->geneticPos[i]  = pos[i];
   }

   *hapDataOut = hapData;
   *mapDataOut = mapData;
}
