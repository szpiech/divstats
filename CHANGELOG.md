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

### Known issues carried forward

- Multi-chromosome input is silently merged into a single coordinate space
  and every row is labelled with the last chromosome seen.
- Unsorted positions are silently accepted and produce wrong windows.
- `--pik k` is intractable for `k >= 3` on realistic data, because tie
  averaging enumerates every `C(t,m)` combination.
- `--map` is required for EHH but the genetic map is never used in any
  computation.
- Output columns are undocumented, and the per-window effective sample size
  is not reported even though it varies with local missingness under the
  default subsampling.

## 1.0.0

Initial version.
