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

## Next, in order

1. **B8/B9** — multi-chromosome input is silently merged into one coordinate
   space and labelled with the last chromosome seen; unsorted positions are
   silently accepted. Both should abort with a diagnostic. Neither changes
   output for well-formed single-chromosome input, so they can land before
   the performance work.
2. **P3** — cache the hypergeometric weight matrix per `(n, H)`. It depends
   only on the pair, not on the window or the data, and the subsampling path
   currently costs about 800x the statistics it feeds (32.9 s against
   0.041 s on 500 windows of 400 haplotypes). This should leave every number
   identical, and it removes the hot path that makes `factln`'s unsynchronised
   memo table worth having.
3. **B3** — closed-form tie averaging in `pi_k2`, replacing exhaustive
   `C(t,m)` enumeration. `--pik k` is currently intractable for k >= 3
   (54 s for 10 windows at 400 haplotypes). The closed form was verified
   against exhaustive enumeration to 9e-13, so values should be identical to
   floating point.
4. **P1** — one-pass VCF reader. Parsing is about 95% of a typical run and
   roughly 4.5x slower than necessary (0.96 s against 0.214 s for a
   one-pass reference reader on the same file).

Then the remaining interface items, of which **U5** (undocumented output
columns) and **U6** (effective per-window sample size not reported) are the
two that affect whether output can be interpreted later. Full list and
measurements in the review report.
