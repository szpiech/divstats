# Branch `divstats-fixes`

Fixes from a code review of `devel` @ `061b966` plus the uncommitted
working-tree changes on top of it. Ordered so that nothing that changes a
reported number lands before there is a suite to measure the change against.

Run `make check` in `src/` after any change on this branch.

## Commits

| commit | what | results change? |
|---|---|---|
| `22e9170` | pre-existing WIP: `--const-n-sub`, macos-arm library path | n/a (baseline) |
| `4d8b090` | **B1** free of uninitialized pointers in `calc_stats` | no |
| `2d53ba6` | golden-file regression suite, `make check` | n/a |
| `d96300c` | **B2** heap over-read in `pi_k2`; leak; uninitialized pointer | no (see below) |
| `904b555` | **B7** `--ehh-part` validated before `--partition` was parsed | no (flag was unusable) |

`22e9170` carries changes that were uncommitted in the working tree when this
branch was cut. To put them back on `devel` as uncommitted work instead:

    git checkout devel && git checkout divstats-fixes -- src/ && git reset

## Why B1 precedes the suite

The suite would normally be the first commit. It cannot be here: a clean
`-O3` build of `22e9170` segfaults on its most basic invocation, so there is
no build to capture goldens from until B1 lands. What establishes that B1
changed nothing is that the binary previously built from this tree — which
does not crash, because the indeterminate pointer values happened to be
benign under the compiler that produced it — passes all goldens captured
from the B1 build.

## Deliberate behavioural differences so far

**None of the commits above changes a computed statistic.** Two are
memory-safety fixes whose effect on output was measured to be nil, and one
enables a flag that previously could not run.

The B2 over-read deserves a precise statement rather than "no change":
`pi_k2` read past the end of `hfs->sortedCount`, and the value it read did
not alias any live key in `count2hap`, so `count()` returned 0 and the
iteration was a no-op. Pre-fix and post-fix builds therefore agree exactly
on every case in the suite. That is a property of this allocator and these
inputs, not a guarantee: had the value aliased a real count, the tie set
would have been built wrongly and `pi_k` silently corrupted. The fix removes
the possibility; it does not correct a number that was previously wrong.

## Suite

21 cases, `tests/run_tests.sh`. Numeric tables are compared column by column
(header exact, integer columns exact, floats to `1e-6` relative), so a
failure names the statistic that moved. Fixtures are derived at run time from
`tests/data/core.vcf.gz`, a 400-site excerpt of `test/test.small.vcf` and the
only input committed — `test/` and `testing/` are untracked, so without it
`make check` cannot run in a clone.

The comparator was negative-controlled before the goldens were trusted:
perturbing a golden float fails and names the column, corrupting a golden
hash fails and prints both, deleting a golden row fails on row count; all
exit non-zero. Every case reproduces itself across three repeated runs, and
`threads-4` is compared against the single-threaded golden so it asserts
thread invariance rather than self-consistency.

Two goldens stand over code with known defects and are expected to move:

- `pik`, `pik-narrow-window`, `pik-multi-k` — `--pik` remains intractable for
  k >= 3 (B3, below). Fixing that replaces exhaustive tie enumeration with a
  closed form; the values should be identical to floating-point, and the
  suite is what will confirm it.
- every case whose input carries missing genotypes — `core.vcf.gz` contains
  `./.` records, so `sites-basic` and friends already exercise the SFS
  subsampler. B4 and B6 below will change those numbers.

## Next, in order

Result-changing fixes. Each needs a golden regeneration with the diff
recorded in the commit message, and collectively they warrant a version bump
since they change published output.

1. **B4** `subsample_sfs` takes the hypergeometric population size as
   `sfs->size`, which is n+1 — the weights behave as if one extra ancestral
   chromosome were in every sample, and the identity projection (H = n) is
   not the identity. Verified against the binary with a six-site fixture:
   divstats reports pi = 2.6, matching the as-coded N = n+1 prediction
   exactly, where the correct projection gives 2.707143. Affects pi, S,
   Tajima's D and Fay & Wu's H on any input with missing data, which is the
   default path. One-character fix, large blast radius.
2. **B6** `segsites()` returns `int`, truncating the fractional subsampled S,
   and that truncated value is what feeds Tajima's D. A `double`-returning
   `s_from_sfs()` already exists in the same file and is never called.
3. **B5** `nCk` materialises the binomial, so it returns `inf` above
   n ~ 1030 on an 8-byte `long double` (arm64 macOS) and the subsampler
   yields nan. Threshold is platform-dependent — roughly 16,000 haplotypes
   on x86-64 Linux — so this reproduces on a cluster and not on a laptop.
   Compute the ratio in log space.
4. **B8/B9** multi-chromosome input is silently merged into one coordinate
   space and labelled with the last chromosome seen; unsorted positions are
   silently accepted. Both should abort with a diagnostic.

Then performance (**P3** cache the hypergeometric weight matrix, **B3**
closed-form tie averaging, **P1** one-pass VCF reader), then the interface
items. Full list and measurements in the review report.

## Correction to the review report

The report claimed the built binary, `.o` files and stale `*.divstats.out`
files were committed to the repo. They were not — `.gitignore` already
covered `*.o`, and the rest were untracked. This branch adds `divstats` and
`*.divstats.out` to `.gitignore` so they stop appearing in `git status`.
