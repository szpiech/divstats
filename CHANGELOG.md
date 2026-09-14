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
- The default `n` is the largest reachable by **99.9% of sites**, so up to 0.1%
  of sites are excluded. Requiring every site to reach it instead lets the
  single worst-covered site set `n` for the whole run: simulating 400
  haplotypes at 2% missing with 0.1% of sites at 50% missing, that rule gives
  `n = 168` against a median site `n` of 392 — a 57% loss of sample for 50
  sites out of 50,000. Startup prints the exclude-nothing value so
  **`--target-n N`** can restore it, or raise `n` further; whenever any site is
  excluded the contributing count appears in **`nSNPsUsed`**. Files with fewer
  than 1000 sites are unaffected, since 0.1% of them rounds to no site.
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

### Changed -- Tajima's D in sub-one-site windows

- **`D` is now undefined when `S < 1`**, written as the na-string. With SFS
  subsampling `S` is an expected count, so `0 < S < 1` is reachable; `D` is
  attenuated there by `sqrt(S) * sqrt(e1/(e1 - e2*(1-S)))`, exactly and
  independently of allele frequency, so a window with a third of an expected
  segregating site reported a `D` shrunk to about 58% of one site's value in
  the same column as fully populated windows. The attenuation is conservative
  — it shrinks toward zero and cannot manufacture a selection signal — but it
  is not a value worth reporting.

  `S`, `pi` and `H` are unchanged and still written for these windows. `D` is
  the only statistic that reads `S`: `H` is the unnormalised *pi* - *theta_H*
  and `pi` does not involve `S`.

  No effect without subsampling, where `S` is a whole number. With
  subsampling this changed 22 of 200 windows at `--winsize 2` and 8 of 80 at
  `--winsize 5` on the 400-site test fixture, and none at `--winsize 100`.

### Added -- usable command line

- **`--help` is grouped, with a usage line and examples.** It was a flat
  alphabetical dump in `std::map` order, so I/O, windowing, statistics and
  missing-data options interleaved, with no usage line and `PREAMBLE` set to
  the empty string. Options now print under *Input*, *Windows*, *Statistics*,
  *Sample size and missing data*, *Output* and *Other*, followed by three
  worked examples and a note on the sample-size and na-string behaviour.

  `param_t` gained `setSectionOrder()` and `setEpilogue()`. The section key is
  the `label` argument `addFlag` already took and which only `printHelp` ever
  read, so no call-site signature changed. A flag whose label is not in the
  declared order still prints, under a trailing *Uncategorised* heading, so
  forgetting to section a new flag cannot hide it. With no order set,
  `printHelp` reproduces the old flat listing.

- **`--help` exits 0.** It exited 1, because `param_t` signals a help request
  with the same `throw 0` it uses for a bad flag and `main` could not tell
  them apart. Both `--help` and `--version` are now handled in `main` before
  parsing, so they work alongside otherwise-invalid arguments.

- **A bare `divstats` prints the usage and exits 1.** It used to reach
  validation and print four errors -- no window mode, window size 0, window
  step 0, no data file -- which is an unhelpful first contact with the tool.

- **Progress on long runs.** Nothing was printed between "Calculating N
  statistics in M windows." and the output file, so a multi-hour run gave no
  sign of life; the progress-bar plumbing in `work_order_t` had been
  commented out. Worker threads now share a counter under a mutex and report
  each 10% of windows to stderr, but only once the run has taken two seconds,
  so short runs stay silent. Measured cost of the lock is +0.3% at one thread
  and +0.4% at four, taking it once per window.

### Changed -- the build

- **The Makefile detects the platform.** It carried three blocks selected by
  commenting and uncommenting, plus a second near-duplicate file
  (`Makefile.lin`) holding the same content with different lines commented,
  so a build change had to be made twice and the default block was wrong for
  most readers. Platform now comes from `uname -s` / `uname -m`; nothing
  needs editing. `Makefile.lin` is deleted.

- **x86 codegen flags are applied only on x86.** `-m64 -mmmx -msse -msse2`
  were unconditional. They do not break an ARM build — clang accepts and
  ignores them — but they emitted "argument unused during compilation" on
  every translation unit.

- **An explicit `-std=c++11`.** There was no `-std=` flag at all, so the
  dialect was whatever the compiler defaulted to, which drifts between
  compiler versions.

- **`make clean`, `install`, `info`.** `clean` removed only `*.o`, leaving
  the binary; it now removes both. `install` takes `PREFIX` and `DESTDIR`.
  `info` prints the detected platform and flags, for bug reports.

- **A missing vendored `libhts.a` says so.** It used to surface as a wall of
  undefined symbols from the linker; `make` now stops with the expected path,
  the detected platform, and the `lib/build_htslib.sh` command that produces
  it.

- **The build is warning-free.** Beyond the x86 flags, seven uses of
  `sprintf` into fixed buffers were deprecated on this toolchain and are now
  `snprintf`. Two of those buffers were genuinely undersized for their
  format: `param_t`'s `char[100]` with `%f` overflows at 1e92 and above
  (101 bytes needed), and `int2str`'s `char[10]` is one byte short of `INT_MIN`
  (`"-2147483648"` plus terminator). Neither is reachable from user input
  today — both take values set in code — but the buffers now fit what the
  types can hold.

- **Stale tracked files removed.** `src/outfile.divstats.out2`, a result file
  in the pre-2.0.0 space-separated header format, and `Makefile.lin`.
  `.gitignore` now also covers `*.divstats.out*` and `*.sweepfinder.out`.

### Added -- run log and output precision

- **Every run writes `<out>.divstats.log`.** Version, timestamp, the full
  command line, input and window parameters, the subsampling target, and a
  partition legend. The legend is the point: statistic columns are suffixed
  `_A`, `_B`, `_C` with nothing recording what each covers, and the answer is
  not recoverable from the table — under `--sites` a partition is a count of
  SNPs, under `--bp` a width in base pairs, and the output does not say which
  mode produced it.

  Written as a sidecar rather than as `#` comment lines above the header, so
  that a bare `pd.read_csv(path, sep='\t')` keeps working. (The code already
  anticipated a log: `work_order_t` carried a commented-out `ofstream *flog`.)

- **`--precision N`.** Statistics were written at the iostream default of 6
  significant digits with no way to change it. The default remains 6, so
  output is unchanged unless asked; 1–17 accepted, and out-of-range is
  refused rather than silently clamped.

### Fixed -- documentation

- **`--partition` units were documented wrongly.** The README said the values
  were percentages under `--sites` and base pairs under `--bp`, summing to
  100 or to `--winsize`. They are in the same units as `--winsize` in both
  modes — SNP counts under `--sites` — and must sum to `--winsize` in both.
  The "100" was an artifact of every example using `--winsize 100`;
  `--sites --winsize 200 --partition 25 50 25` is refused with "Window
  partitions sum to 100 but must sum to 200 instead."

### Changed -- parsing

- **The genotype matrix is accumulated site-major and transposed once.** It
  was appended per haplotype — one `push_back` into each of `nhaps` separate
  buffers at every record, so 20 million scattered appends on a
  400-haplotype × 50,000-site file, each landing on a different cache line,
  with `nhaps` independent growth schedules. Records now write one contiguous
  run of `nhaps` bytes (one `resize` per record, then plain stores), and the
  matrix is transposed at the end in cache-sized blocks of loci. divstats'
  own share of the parse fell from 0.133 s to 0.039 s, a factor of 3.4.

- **`--threads` is handed to htslib's decoder.** This parallelises BGZF block
  inflation, so it does something for a bgzipped VCF or a BCF and nothing for
  plain gzip. Measured at 3–7% of the parse: the dominant cost in the text
  path is `vcf_parse`, which is serial.

Parse time, 400 haplotypes × 50,000 sites, 4 threads: plain gzip VCF 0.347 →
0.288 s, bgzipped VCF 0.331 → 0.272 s, BCF 0.230 → 0.148 s. End to end over
500 windows: 13.5%, 15.5% and 25.4% faster respectively. Output is
byte-identical in all three formats at one and four threads, and under
`--hemi`.

The largest lever is the input format, not the code: BCF cuts htslib's share
from 0.183 s to 0.032 s, because it removes the text parse entirely. This is
now documented in the README.

### Known issues carried forward



## 1.0.0

Initial version.
