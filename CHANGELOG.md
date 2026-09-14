# Changelog

## 2.0.0

**Statistics computed on data with missing genotypes have changed value.**
Do not pool output from 1.0.0 with output from 2.0.0. If you have published
or archived divstats results, they were computed with the defects below and
should be regenerated.

The affected statistics are the SFS-derived ones — `pi`, `S`, Tajima's `D`,
Fay & Wu's `H`, and their partition columns. They change only when the SFS
subsampler runs, which is whenever an input has missing genotypes and
`--no-sfs-sub` was not given (subsampling is on by default). EHH statistics
(`ehh*`) and `pik` are unaffected. Output on data with no missing genotypes
at all is unchanged.

Measured on the regression fixtures, comparing each fix against the commit
before it, as maximum relative deviation per column:

| fix | π | S | D | H |
|---|---|---|---|---|
| hypergeometric population size | 18.8% | 16.7% | 69.6% | 6792% |
| fractional segregating sites | — | 10.0% | 58.4% | — |
| log-space subsampling weights | 12.6% | 14.1% | 55.4% | 534% |

The very large `D` and `H` figures are windows where the statistic sits near
zero, so a small absolute shift is a large relative one. Nucleotide diversity
moving by up to roughly 19% is the more representative headline.

### Fixed — correctness

- **Hypergeometric population size was off by one.** `subsample_sfs` took the
  population size from `sfs->size`, but the array spans indices `0..n`, so
  that value is `n+1`. Every weight was computed as if one extra ancestral
  chromosome were present, and projecting `n` onto `H = n` was not the
  identity. Fixing it requires two coupled changes — the population size and
  the loop bound, which was expressed in terms of it.
- **Fractional segregating-site counts were truncated to `int`.** When the
  SFS is subsampled its bin counts are fractional expectations, so `S` is
  fractional. The truncation discarded up to a whole site from the reported
  `S` column and from the value fed into Tajima's `D` as both numerator and
  variance term.
- **Subsampling weights are now computed in log space.** Forming each
  binomial coefficient separately was wrong twice over. It overflowed to
  `inf` for large cohorts — the *product* of two mid-range coefficients
  overflows long before `C(n,H)` does, and `inf/inf` is `nan`, so every SFS
  statistic came back `nan` at 1200 haplotypes with 5% missing data. It also
  invented probability mass: `nCk` rounds, and outside the hypergeometric
  support it returns 1 where the true coefficient is 0, so one 100-site
  window projected to 101.65 sites. The overflow threshold is
  platform-dependent (`long double` is 8 bytes on arm64 macOS, 80-bit on
  x86-64 Linux), which is why this reproduced on some machines and not
  others.

### Fixed — crashes and memory safety

- **A clean build of 1.0.0 segfaults on its most basic invocation.**
  `calc_stats` declared `sfs`, `partition_sfs` and `pik_hfs` uninitialized,
  assigned them inside conditional branches, and released all three
  unconditionally on every window iteration. The binary distributed with
  1.0.0 does not crash only because the indeterminate values happened to be
  benign under the compiler that produced it.
- **Heap over-read in `pi_k2`.** The loop collecting haplotype frequency
  classes was bounded by `k`, but `hfs->sortedCount` holds `hfs->size`
  entries — the number of distinct *counts*, which is generally far smaller
  than the number of distinct haplotypes that `k` is checked against. Any
  window where the distinct-count cardinality fell below `k` read past the
  allocation. Also fixes a `string[]` leak on one early-return path and a
  test of an uninitialized pointer in the same function.
- **The Makefile declared no header dependencies**, so editing a header left
  stale objects linked against the old declarations. This is not cosmetic:
  changing `segsites()` from `int` to `double` against a stale object
  produced a binary that read the return value out of the integer register
  and printed `S = 5.66e+07`. If you have built this tree incrementally
  after a header edit, assume your binary was mislinked.

### Fixed — interface

- **`--ehh-part` could never be used.** Its consistency check ran before the
  block that parses `--partition`, so it tested a variable that still held
  its initializer and always failed. The error message also named the wrong
  flag.
- The `--partition` count error reported the wrong variable in its "maximum
  allowed" message.

### Added

- **`make check`** — a golden-file regression suite (`tests/run_tests.sh`),
  23 cases, comparing numeric output column by column so a failure names the
  statistic that moved. Fixtures derive at run time from three small
  committed inputs, so it runs in a fresh clone; `test/` and `testing/` are
  untracked and absent from clones.

### Added -- input formats and performance

- **BCF support**, plus a single-pass reader. VCF, bgzipped VCF and BCF now
  read through htslib, which detects the format from the file itself. The old
  pair of functions read every input four times -- each counted records, closed
  the file and reopened it, and the map pass re-parsed CHROM and POS off lines
  the genotype pass had already read. Parsing is about 2.6x faster
  (1.02 s -> 0.40 s on 50,000 sites x 400 haplotypes) with byte-identical
  output. `GT` may now sit anywhere in `FORMAT`; it was previously assumed to
  be first, and `FORMAT=DP:GT` failed with a misleading allele-coding error.
- **Multi-chromosome and unsorted input are refused**, naming the offending
  record, instead of being analysed incorrectly.
- **Subsampling is about 44x cheaper.** The projection evaluated a weight for
  every input frequency bin, almost all of which are empty; it now visits only
  the bins that hold sites, which makes the cost independent of sample size.
  On 500 windows of 400 haplotypes with 2% missing data, the subsampling cost
  above the no-subsampling baseline falls from 12.23 s to 0.28 s -- total
  runtime 13.31 s to 1.36 s. Output is bit-identical.
- **`--pik k` is now usable for k >= 3.** Tie averaging enumerated every
  C(t,m) choice of tied haplotypes, which is intractable when most haplotypes
  in a window are unique. It has a closed form, so the cost no longer depends
  on k: on 40 haplotypes over 10 windows, k=5 goes from 10.09 s to 0.063 s,
  with identical output.

### Changed -- build

- **GSL is no longer required.** `gsl_combination` was its only use, in the
  tie enumeration the closed form replaces. 52 MB of vendored static libraries
  and 1.4 MB of headers are gone from the repository.
- **htslib is vendored** in its place (4.5 MB), configured without
  bzip2/lzma/libcurl/libdeflate so it needs nothing beyond the zlib divstats
  already linked. Only `lib/macos-arm/libhts.a` is committed;
  `lib/build_htslib.sh` reproduces it for other platforms.
- Malformed input now exits 1 with its diagnostic. The readers signal failure
  by throwing, but only command-line parsing was wrapped, so these aborted on
  SIGABRT (exit 134) after printing the error.

### Added -- genetic distance

- **`--ehh-cm w [w…]`** places EHH sub-windows by genetic distance, in the
  map's own units. A fixed genetic width covers fewer base pairs where
  recombination is high, so values are comparable between regions in a way a
  width in base pairs is not. Columns are `ehhcm_<w>`. Not computed per
  partition.
- **`--ehh` no longer requires a genetic map.** It demanded `--map` or
  `--pmap` for any EHH calculation and then placed its sub-windows by physical
  position regardless, so the map was parsed, stored and never read -- users
  produced a recombination map for a computation that ignored it. Physical
  placement is now the default and `--ehh-cm` is the only consumer of
  `geneticPos`. `--pmap` is kept but is only needed to force physical
  placement while a map is loaded for something else, and is refused together
  with `--ehh-cm`.

### Changed -- one sample size for the whole run

- **Windows are now comparable.** Each window used to be projected to its own
  minimum observed sample size, so `n` varied with local missingness while the
  output gave no way to tell. Every window is now projected to a single `n`,
  and that `n` is reported in a new **`nhaps`** column. This changes π, S, *D*
  and *H* for any input with missing genotypes -- on the test fixture,
  74 of 200 windows move at `--winsize 2`, by up to a factor of 2.5 where the
  old per-window `n` happened to be high.
- The default `n` is the largest every site can reach, which excludes nothing.
  Because one badly-covered site sets it, startup reports what raising it would
  cost in sites, and **`--target-n N`** projects to `N` instead, excluding
  sites too sparse to reach it and reporting the contributing count in
  **`nSNPsUsed`**.
- **`--window-n-sub`** restores the previous per-window behaviour.
  **`--const-n-sub` is now a no-op**: it selected what is now the default, and
  its regression case shares the default's expected output to assert that.

### Added -- `--version`

- `divstats --version` prints the version to stdout and exits 0. It was
  previously reachable only in the stderr banner, mixed with progress output,
  which conda recipes and pipeline version-capture cannot use.

### Fixed -- output and parsing

- **Undefined statistics are `nan`, not `-999`.** A numeric sentinel is
  absorbed by anything that averages a column. `--na-string` chooses the token
  (`NA`, empty, or `-999` to reproduce 1.x).
- **The header is tab-separated** like the data rows. It joined statistic names
  with spaces, so the header had 6 tab-delimited fields where rows had 9 and
  both `read.table(header=TRUE)` and `pandas.read_csv(sep='\t')` mis-aligned.
- **The header no longer depends on a worker thread.** Column names were
  appended while thread 0 processed window 0, so a run with zero windows wrote
  a header with no statistic columns.
- **SNPs sitting exactly on a window boundary are included.** With `--bp
  --winsize 1000` and SNPs at 999, 1999, ... every window previously reported
  zero SNPs.
- **Tajima's D no longer divides 0 by 0** where a window has no segregating
  sites.
- **`--sweepfinder` reports the sample size observed at each site**, not the
  full cohort size, and writes `<out>.sweepfinder.out` rather than stdout. It
  is now documented in `--help`.
- **Command-line parsing**: boolean flags are set rather than toggled; `""`,
  `"-"` and `"."` are no longer accepted as numbers; integer overflow is
  detected; and every bad flag is reported in one run instead of only the
  first.

### Known issues carried forward


- The partition columns are labelled `_A`, `_B`, `_C` with no record anywhere
  of which base-pair span each letter covers, so a partitioned output file
  cannot be decoded without the original command line.

## 1.0.0

Initial version.
