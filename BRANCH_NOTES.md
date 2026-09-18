# Branch `divstats-fixes`

Fixes from a code review of `devel` @ `061b966` plus the uncommitted
working-tree changes that sat on top of it. Ordered so that nothing which
changes a reported number lands before there is a suite to measure the
change against.

Run `make check` in `src/` after any change on this branch.

## Commits

| commit | what | changes reported numbers? |
|---|---|---|
| `22e9170` | pre-existing WIP: `--const-n-sub`, macos-arm library path | n/a (baseline) |
| `4d8b090` | **B1** free of uninitialized pointers in `calc_stats` | no |
| `2d53ba6` | golden-file regression suite, `make check` | n/a |
| `d96300c` | **B2** heap over-read in `pi_k2`; leak; uninitialized pointer | no |
| `904b555` | **B7** `--ehh-part` validated before `--partition` was parsed | no (flag was unusable) |
| `031636b` | these notes | n/a |
| `29a3ba0` | **B4** hypergeometric population size was n+1 | **yes** |
| `fb68b0f` | **B6** fractional S truncated to int, and fed to Tajima's D | **yes** |
| `cbc75df` | harness: nan counted as equal; cleanup trap deleted the repo | n/a |
| `3a4d802` | **B5** subsampling weights in log space | **yes** |
| `ed592ce` | these notes, rewritten with the measured deviations | n/a |
| `100c335` | correct the overflow-threshold note in the large-n comment | n/a |
| `8fb6226` | complete the commit table in these notes | n/a |
| `9e61816` | bump to 2.0.0, add CHANGELOG | n/a |
| `7f2c81e` | **P3** project the SFS from nonzero bins only, 44x | no (bit-identical) |
| `56a45d1` | **B3** closed-form tie averaging; drop GSL | no (bit-identical) |
| `4ece3b6` | **BCF** + one-pass htslib reader; **B8 B9 B11 P1** | no (byte-identical) |
| `9700233` | README, CHANGELOG, BRANCH_NOTES for 2.0.0 | n/a |
| `8a48b0a` | **B13 B14** header built in main, tab-separated | header line only |
| `671deae` | **B15 B16** undefined statistics are nan; guard D | `-999` -> `nan`, 69 values |
| `1a6728a` | **B12** SNPs exactly on a window boundary | no golden moved |
| `3f48fa7` | **B10 B19 U4** sweepfinder n and destination; parser | sweepfinder golden |
| `6866ffe` | **B17** `--ehh-cm` genetic placement; `--ehh` drops the map | no (new columns only) |
| `0183804` | **B18** `releaseHapData` never freed the struct | no |
| `36bd320` | **U3 U6** `--version`; one sample size per run, `nhaps` column | yes, with missing data |
| `ae582f4` | default target = n reachable by 99.9% of sites | no (fixtures <1000 sites) |
| `7f8bc2a` | **P8** flags read once, not per window | no (byte-identical) |
| `fa9f913` | **P4** share EHH sub-window spectra with `--ehhk` | no (byte-identical) |
| `ba0baea` | **P5** one copy per haplotype string | no (byte-identical) |
| `c942abf` | **P9** count uniques in the existing pass | no (byte-identical) |
| `bba3052` | **P10** contiguous results table, no per-row flush | no (byte-identical) |
| `c4d0251` | Tajima's D undefined below one expected segregating site | yes, D only |
| `573c70d` | **U1 U2 U7** grouped `--help`, usage on bare run, progress | no (byte-identical) |
| `48f8b15` | **U10** platform-detecting Makefile; warning-free build | no (byte-identical) |
| `a5de078` | **U5 U9** `<out>.divstats.log` sidecar; `--precision` | no (byte-identical) |
| `7504fdd` | **P2** site-major genotype fill + blocked transpose; htslib thread pool | no (byte-identical) |
| `036ae89` | **P12** `--bp` windows start at the data, not coordinate 0 | yes: empty leading rows gone |
| `a3d5560` | **P11** alleles stored two bits per site | no (byte-identical) |
| `37b67d9` | **P6 P7** hoisted spectrum lookups, const-ref strings and maps | no (byte-identical) |
| `48034cc` | MSYS2/MinGW64 as a platform: explicit uname match, `$(EXE)`, mingw64 archive | no |
| `081a647` | LF output on every platform (binary-mode streams); suite finds `divstats.exe` | Windows only |
| `d09e414` | MSYS2 CI job; `.gitattributes` pins the tree to LF | no |

For the count and anything added after the rows above, ask git rather than
trusting this table:

    git log --oneline devel..HEAD
    git rev-list --count devel..HEAD

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

## What changed in the numbers

B1, B2 and B7 change nothing: two are memory-safety fixes whose effect on
output was measured to be nil, and one enables a flag that could not run.
The B2 over-read deserves a precise statement rather than "no change":
`pi_k2` read past the end of `hfs->sortedCount`, and the value it read did
not alias any live key in `count2hap`, so `count()` returned 0 and the
iteration was a no-op. That is a property of this allocator and these
inputs, not a guarantee — had the value aliased a real count, the tie set
would have been built wrongly and `pi_k` silently corrupted.

**B4, B6 and B5 all change published output.** Together they warrant a
version bump and a note in the release history. Anyone holding divstats
output computed with subsampling on (the default) and any missing genotypes
should expect different numbers. Deviations were measured by running the
preceding commit's binary and the new one against identical inputs:

| fix | cases moved | largest per-column deviation |
|---|---|---|
| B4 | 14 of 22 | pi 18.8%, S 16.7%, D 69.6%, H 6792%; partition D_B 960% |
| B6 | 12 of 22 | S 10.0%, D 58.4%; partition S_C 29.3%, D_A 63.0%, D_C 33214% |
| B5 | 14 of 23 | pi 12.6%, S 14.1%, D 55.4%, H 534%; partition D_A 30.5% |

The very large H and D figures are windows where the statistic sits near
zero, so a small absolute shift is a large relative one. The honest headline
is that nucleotide diversity moves by up to roughly 19%, 13% and 13% at the
three steps.

In every case the statistics that move are exactly the SFS-derived ones.
`missing-nosub` (subsampling disabled), the three pure-EHH cases and the
three `pik` cases are byte-identical across all three fixes, and no `ehh_*`
column changes anywhere. The `ehh-part` case does move under B4 and B5,
because it reports `pi` columns alongside its `ehh` columns.

## Two corrections to the review report

1. The report gave the B4 fix as the one-character `int n = sfs->size - 1;`.
   That is wrong and would introduce silent data loss: the loop bound is
   coupled to `n`, and leaving `i < n - 1` drops the top segregating class.
   On a fixture with a site at `c == n-1` the site vanishes entirely
   (projected mass 3.0 from 4 input sites). The fix is two lines.
2. The report treated B5 as purely an overflow guard. It is also a
   correctness fix: `nCk` rounds, and outside the hypergeometric support it
   returns 1 instead of 0 for `(a,b)` in `{(0,1),(0,2),(1,2)}`, inventing
   probability mass. One 100-site window projected to 101.65 sites; another
   to 103.30. This affects far smaller cohorts than the overflow does.

The report also claimed the built binary, `.o` files and stale
`*.divstats.out` files were committed. They were not — `.gitignore` already
covered `*.o` and the rest were untracked.

## Suite

23 cases, `tests/run_tests.sh`. Numeric tables are compared column by column
(header exact, integer columns exact, floats to `1e-6` relative), so a
failure names the statistic that moved. Fixtures derive at run time from
three committed inputs: `tests/data/core.vcf.gz` (a 400-site excerpt of
`test/test.small.vcf`), `subsample-toy.vcf.gz` (hand-computable, 6 sites) and
`large-n.vcf.gz` (1200 haplotypes, guards the overflow). Committing small
fixtures is deliberate: `test/` and `testing/` are untracked, so without them
`make check` cannot run in a clone.

Four negative controls, all confirmed to fail before the goldens were
trusted: a perturbed golden float (names the column), a corrupted golden
hash, a deleted golden row (row count), and `nan` against a number. The
fourth exists because the first three all passed while the comparator was
treating `nan` as equal to every number — a negative control only tests the
failure mode it encodes. `threads-4` is compared against the
single-threaded golden, so it asserts thread invariance rather than
self-consistency.

## Incident, 2026-09-12

The suite's cleanup trap ran `rm -rf` on the repository root. `mktemp -d`
failed because `TMPDIR` was set to a relative, nonexistent path, leaving
`SCRATCH` empty; a line added to absolutise it,
`SCRATCH="$(cd "$SCRATCH" && pwd -P)"`, then resolved it to the current
directory, because bash treats `cd ""` as a successful no-op. The trap's
unguarded `rm -rf "$SCRATCH"` deleted the working tree, most of `.git`, and
the untracked `test/` and `testing/` directories.

Recovered from Dropbox file history plus a `git archive` snapshot: all
commits intact, every tracked file byte-identical, and `test/` and `testing/`
restored. Confirmed by re-deriving `tests/data/core.vcf.gz` from the restored
`test/test.small.vcf` — the sha256 matches the committed fixture — and by a
clean rebuild passing the suite. `git fsck` is clean.

Guarded in `cbc75df`: `mktemp` failure is now fatal before anything is
registered for deletion, the resolved path must match
`*/divstats-tests.??????`, and cleanup re-checks that shape and refuses any
other path. If you run the suite and it exits 2 complaining about the scratch
directory, that guard is doing its job — check `TMPDIR`.

## Performance, measured

Each change was checked for output equivalence before its speed was recorded.

> **Benchmark caveat.** The numbers for P1, P3, P5 and P8 below were taken on
> a `bench/big400miss.vcf.gz` whose header never declared `##FORMAT=<ID=GT>`,
> so htslib fell back to a string-typed GT on every record. The fixture was
> regenerated with the declaration during P2; those rows are internally
> consistent but are not directly comparable with the P2 figures.

| change | effect | equivalence |
|---|---|---|
| P3 sparse SFS projection | subsampling cost above baseline 12.23 s -> 0.28 s (44x); total 13.31 s -> 1.36 s (9.8x) | bit-identical; suite passes at `--rtol 1e-15`, `cmp` clean over 500 windows |
| B3 closed-form tie averaging | `--pik 5` on 40 haplotypes, 10 windows: 10.09 s -> 0.063 s; flat in k | identical output at every k tested |
| P1 one-pass htslib reader | parse-dominated run 1.02 s -> 0.40 s (2.6x) | byte-identical over 500 windows; all 25 prior goldens unchanged |

P3's implementation deliberately differs from what the review proposed. The
review called for caching the weight matrix per `(n,H)`; `H` is the
per-window minimum sample size, so it varies window to window and the key
space is the product of observed `n` and `H` values. At 1200 haplotypes each
matrix is ~11 MB, so a few hundred live pairs would reach gigabytes. Skipping
empty input bins gives a comparable speedup with no memory growth.

## A coverage gap worth remembering

The three `pik` cases all reached the tie-averaging branch only with
`t=2, m=1` -- one haplotype drawn from a class of two -- which was
established by instrumenting the branch with a counter, not assumed. In that
trivial case the closed form's third term has coefficient zero, so the
byte-identical results the suite reported said nothing about it.
`pik-tie-m2` and `pik-tie-m3` were added to reach `t=4, m=2` and `t=4, m=3`.
The general lesson matches the one from the `nan` comparator hole: a passing
test says nothing about a path it does not execute, and the cheapest way to
find out which paths a case reaches is to instrument and look.

## Next, in order

Every bug in the review is fixed (19 of 19), as are U3, U4, U6, U8, U11 and
the P1/P3/P4/P5/P8/P9/P10 performance items. What remains:

1. **P2, P6, P7** -- `hamming_dist_str` takes both strings by value and is
   called from the O(k^2) inner loops; the map lookups in those loops are
   O(log n) tree walks with O(L) string comparisons per pair.
2. **U1, U2, U7** -- usage on bare invocation, grouped `--help`, progress
   reporting on long runs. `--help` is currently an alphabetical dump in
   `std::map` order with no usage line and no examples.
3. **U9, U10** -- output precision flag; a Makefile keyed on `uname`
   rather than commented blocks; and the stale `.o` files and result files
   still checked in.
4. **P11** -- 2-bit genotype packing, only if cohort-scale runs are a target.

Full list and measurements in the review report.

## What the P8/P4/P5/P9/P10 tier actually bought

Measured on a 400-haplotype x 50,000-site file with 2% missing genotypes, best
of three to seven runs, against the build immediately before the tier. Every
row is byte-identical output.

| item | configuration | before | after | |
|---|---|---|---|---|
| P8 | 25,000 windows, 1 thread | 0.762 s | 0.769 s | noise |
| P8 | 25,000 windows, 4 threads | 0.578 s | 0.533 s | 7.8% |
| P4 | `--ehh` **and** `--ehhk` | 0.497 s | 0.420 s | 15.5% |
| P4 | `--ehh` alone | 0.416 s | 0.418 s | unchanged |
| P5 | 10,000-window slide, `--ehh 100` | 1.12 s | 0.57 s | 49.4% |
| P5 | 10,000-window slide, `--ehh 20 50` | 2.07 s | 1.60 s | 22.7% |
| P9 | same, on top of P5 | 22.7% | 23.3% | noise |
| P10 | 9,982 rows | 6.221 s | 6.066 s | 2.5% |

Three of the five are not really performance changes and the notes should say
so rather than let the table imply otherwise:

- **P8** is a correctness fix. 200,000 map lookups are nothing next to the
  arithmetic they guarded, which is why the single-thread column is flat. What
  it removes is every worker thread calling a non-const member of one shared
  `param_t` for the whole run -- undefined behaviour, not merely slow. The
  4-thread gain is the contention disappearing.
- **P9** contributes one to two points over P5 alone, inside run-to-run
  variation. It stands on removing a per-call red-black tree from the
  per-window path and on a latent out-of-bounds read of `array[0]` before the
  size was checked.
- **P10** is 0.8-2.5%, as a low-severity finding should be. It also supplies
  the `delete` the results table never had -- allocated once per run and held
  to exit, so never a growing leak, but 160 MB for a million windows with 20
  statistics.

**P5 is the one with real reach**, and measuring it took three attempts. On a
coarse scan the spectra cost 0.015 s inside a 0.33 s run that is otherwise
parsing, and the two builds differed by noise -- one configuration looked 3%
*slower*. Subtracting a "parse floor" made it worse, yielding negative phase
times, because a floor measured at a different window size is not a floor. The
effect only appears when `hfs_window` is called enough to matter, which for
EHH means a sliding scan. The lesson is the one P8 repeats from the other
direction: decide what configuration exercises the code you changed before
reading any number off it.

## Two places where verification changed the answer

Worth keeping because in both cases the first result was wrong in the
reassuring direction.

**A guard that cannot fire.** The radicand check added next to Tajima's D
looked like it was catching a real nan. Computing `e1` and `e2` over
n = 4..5000 showed `e1 > e2` throughout, so the radicand is positive for any
S; at n = 2 and 3 both coefficients are exactly 0, but pi == S/a1 identically
there, so the answer is nan either way. Instrumenting the branch over ~15,000
fractional-S evaluations reached it zero times. The guard stays, but the
comment now says it is defensive rather than claiming a fix -- the first
version of that comment asserted the opposite and was wrong.

**A sanitizer that was not running.** The `releaseHapData` leak fix appeared
to be confirmed by AddressSanitizer reporting no leak. It was not: LSan does
not run on macOS arm64, and a control program with a deliberate leak produced
no report either. The fix rests on inspection and on matching the other three
release functions. Re-check on Linux with `detect_leaks=1`.

**Windows is unproven until CI runs it.** There is no Windows machine in the
environment the port was written in, so every claim about the mingw64 leg is
either platform-independent or verified by simulating `uname`: all eight
`L_PATH` branches by overriding `UNAME_S`/`UNAME_M` on the make command line,
the build script's guesser through a fake `uname`, `check_htslib.sh` against a
non-ELF archive and an unreadable one, and the `.exe` fallback by pointing
`--bin` at a suffixed copy. What has never executed is a MinGW compile, a
native-Windows run, or the CRLF behaviour the binary-mode streams exist to
prevent. The `windows` job in CI is the test; read it before trusting the
platform.

Two things were found while writing it that were defects in their own right,
both noted in the commits: the Windows leg of the Makefile was the catch-all
`else`, so every unrecognised platform got Windows library paths; and its
`OPT ?=` could never fire, because `OPT ?= -O3` is resolved above the platform
block -- so `-static-libgcc` and `-DPTW32_STATIC_LIB` had been dead since they
were written.
