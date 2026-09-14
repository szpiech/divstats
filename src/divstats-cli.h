#ifndef __DIVSTATS_CLI_H__
#define __DIVSTATS_CLI_H__

//2.0.0 rather than 1.1.0: pi, S, Tajima's D and Fay & Wu's H all changed
//value on data with missing genotypes (see CHANGELOG.md). For a tool whose
//interface IS its numbers, that is a breaking change, and a major version is
//the only signal strong enough to stop someone pooling 1.0.0 and 2.0.0
//output in one analysis.
const string VERSION = "2.0.0";

const string PREAMBLE = "\n\
divstats v" + VERSION + " -- window-based diversity statistics\n\
\n\
USAGE:  divstats --vcf <file> (--sites | --bp) --winsize <n> --winstep <n> \\\n\
                 <statistic> [<statistic> ...] --out <prefix>\n\
\n\
Reads VCF, bgzipped VCF, BCF or TPED and writes one row per window to\n\
<prefix>.divstats.out. Every statistic is optional; ask for at least one.\n\
Input must be a single chromosome, sorted by position.";

const string EPILOGUE = "\n\
EXAMPLES\n\
--------\n\
\n\
  Nucleotide diversity, segregating sites, Tajima's D and Fay & Wu's H in\n\
  non-overlapping 100-SNP windows:\n\
\n\
    divstats --vcf chr2.vcf.gz --sites --winsize 100 --winstep 100 \\\n\
             --pi --s --d --h --out chr2\n\
\n\
  The same in 10 kb windows sliding by 2 kb:\n\
\n\
    divstats --vcf chr2.vcf.gz --bp --winsize 10000 --winstep 2000 \\\n\
             --pi --s --d --h --out chr2\n\
\n\
  EHH in 20- and 50-SNP sub-windows, on 8 threads:\n\
\n\
    divstats --vcf chr2.bcf --sites --winsize 200 --winstep 200 \\\n\
             --ehh 20 50 --threads 8 --out chr2\n\
\n\
NOTES\n\
-----\n\
\n\
  With missing genotypes every window is projected to one sample size, fixed\n\
  for the whole run and reported in the nhaps column; see --target-n. Windows\n\
  where a statistic is undefined are written as nan; see --na-string.";

const string ARG_THREADS = "--threads";
const int DEFAULT_THREADS = 1;
const string HELP_THREADS = "The number of threads to spawn during computations.";

// I/O flags
const string ARG_FILENAME_TPED = "--tped";
const string DEFAULT_FILENAME_TPED = "__hapfile1";
const string HELP_FILENAME_TPED = "A TPED file containing haplotype and map data.\n\
\tVariants should be coded 0/1";

const string ARG_FILENAME_POP1_VCF = "--vcf";
const string DEFAULT_FILENAME_POP1_VCF = "__vcffile1";
const string HELP_FILENAME_POP1_VCF = "A VCF file containing haplotype data.\n\
\tVariants should be coded 0/1";

const string ARG_FILENAME_MAP = "--map";
const string DEFAULT_FILENAME_MAP = "__mapfile";
const string HELP_FILENAME_MAP = "A mapfile with one row per variant site.\n\
\tFormatted <chr#> <locusID> <genetic pos> <physical pos>.";

const string ARG_HEMI = "--hemi";
const bool DEFAULT_HEMI = false;
const string HELP_HEMI = "Input data is hemizygous (e.g. chrY).";

const string ARG_TARGET_N = "--target-n";
const int DEFAULT_TARGET_N = 0;
const string HELP_TARGET_N = "Project every window's spectrum to exactly this many haplotypes.\n\
\tSites observed in fewer than this are excluded (they cannot be\n\
\tprojected upward), and the count that contributed is reported in the\n\
\tnSNPsUsed column. The default is the largest n every site in the file\n\
\tcan reach, which excludes nothing but is set by the single worst\n\
\tcovered site; the startup message reports what raising it would cost.";

const string ARG_WINDOW_N_SUB = "--window-n-sub";
const bool DEFAULT_WINDOW_N_SUB = false;
const string HELP_WINDOW_N_SUB = "Project each window to ITS OWN minimum sample size rather than to one\n\
\tvalue shared by the whole run. Maximises each window's n in isolation,\n\
\tbut n then varies with local missingness and pi, S, D and H are NOT\n\
\tcomparable between windows. This was the default before 2.0.0.";

const string ARG_PRECISION = "--precision";
const int DEFAULT_PRECISION = 6;
const string HELP_PRECISION = "Significant digits written for each statistic. The default, 6, is\n\
\tthe iostream default and what every previous version emitted.";

const string ARG_VERSION = "--version";
const bool DEFAULT_VERSION = false;
const string HELP_VERSION = "Print the version to stdout and exit.";

const string ARG_NA_STRING = "--na-string";
const string DEFAULT_NA_STRING = "nan";
const string HELP_NA_STRING = "Token written for statistics that are undefined in a window\n\
\t(too few sites, no segregating sites, fewer haplotypes than --pik k).\n\
\tPrevious versions wrote -999, a numeric value that any downstream\n\
\tmean() or quantile() absorbs as data.";

const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "The basename for all output files.";

const string ARG_2_SWEEPFINDER = "--sweepfinder";
const bool DEFAULT_2_SWEEPFINDER = false;
const string HELP_2_SWEEPFINDER = "Write SweepFinder2 input to <out>.sweepfinder.out and exit.\n\
\tColumns are position, derived count, haplotypes observed at that site,\n\
\tand folded. folded is always 0, i.e. the ALT allele is assumed to be\n\
\tthe derived one -- divstats has no outgroup and cannot verify this.\n\
\tSites with no called genotype are omitted.";

// Window control flags
const string ARG_BP = "--bp";
const bool DEFAULT_BP = false;
const string HELP_BP = "Use bps for window sizes.";

const string ARG_SITES = "--sites";
const bool DEFAULT_SITES = false;
const string HELP_SITES = "Use sites for window sizes.";

const string ARG_WINSIZE = "--winsize";
const int DEFAULT_WINSIZE = 0;
const string HELP_WINSIZE = "The window size within which to calculate diversity statistics.";

const string ARG_WINSTEP = "--winstep";
const int DEFAULT_WINSTEP = 0;
const string HELP_WINSTEP = "The sliding window step size.";

const string ARG_PARTITION = "--partition";
const int DEFAULT_PARTITION = 0;
const string HELP_PARTITION = "Partition the sliding window into non-overlapping sub windows\n\
\tof varying sizes within which all statistics (except EHH-based ones) are calculated separately.\n\
\te.g. For a sliding window of 100kb, --partition 25000 50000 25000 would instruct\n\
\tdivstats to calculate statistics separately within a central 50kb and within the\n\
\ttwo flanking 25kb regions for each 100kb window. Partitions must add up to --winsize.\n\
\tSet to 0 to simply calculate within the entire window.";
const int MAX_PARTITION = 20;

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
const string HELP_EHH = "A list of window sizes within which to calculate EHH.\n\
\tThese subwindows will be centered on the middle of the current\n\
\twindow and may not be larger than --winsize.\n";

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

// Other flags
const string ARG_EHH_PART = "--ehh-part";
const bool DEFAULT_EHH_PART = false;
const string HELP_EHH_PART = "Calculates EHH/EHHK in any partitions of the main window.\n\
\tTo be distinguished from --ehh, which calculates EHH in sub-windows\n\
\tcentered on the main window.  Requires --ehh/--ehhk to be set.";

const string ARG_NO_SFS_SUB = "--no-sfs-sub";
const bool DEFAULT_NO_SFS_SUB = false;
const string HELP_NO_SFS_SUB = "Do not subsample the SFS to handle missing data.\n\
\tEffectively treats missing data as 0/0 in the SFS. Can speed up computation of\n\
sfs-based statistics substantially.";

const string ARG_CONST_N_SUB = "--const-n-sub";
const bool DEFAULT_CONST_N_SUB = false;
const string HELP_CONST_N_SUB = "Accepted for compatibility and does nothing: a constant sample size\n\
\tacross windows is the default from 2.0.0. Use --window-n-sub for the\n\
\tprevious per-window behaviour, or --target-n to choose the value.";

const string ARG_PMAP = "--pmap";
const bool DEFAULT_PMAP = false;
const string HELP_PMAP = "Place EHH subwindows by physical distance. This is now the default\n\
\twhen no --map is given, so the flag is only needed to force physical\n\
\tplacement while a map is loaded for something else.";

const string ARG_EHH_CM = "--ehh-cm";
const double DEFAULT_EHH_CM = 0;
const string HELP_EHH_CM = "A list of subwindow widths in GENETIC distance, in whatever units the\n\
\t--map file uses (normally cM), within which to calculate EHH. Requires\n\
\t--map. Subwindows are centred on the genetic midpoint of the current\n\
\twindow. Unlike --ehh, whose values are SNP counts under --sites and\n\
\tbase pairs under --bp, these widths are comparable between regions of\n\
\tdiffering recombination rate. Decimals are expected: 1 cM is a very\n\
\tlarge window.";


#define NOPTS 8

//Column order follows this array, not the order flags appear on the command
//line. --ehh-cm sits next to the other EHH variants.
const string STATS[NOPTS] = {ARG_PI, ARG_PIK, ARG_SEGSITES, ARG_EHH, ARG_EHHK, ARG_EHH_CM, ARG_TAJ_D, ARG_FAY_WU_H};

#endif