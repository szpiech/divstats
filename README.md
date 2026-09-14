# divstats

Window-based diversity statistics for genetic data analysis.

Computes nucleotide diversity, segregating sites, Tajima's *D*, Fay & Wu's
*H*, extended haplotype homozygosity and *k*-haplotype variants of π and EHH,
in windows defined either by a number of SNPs or by physical distance.

## Version

    divstats --version        # prints e.g. 2.0.0 on stdout, exits 0

## Build

    cd src && make

Requires a C++ compiler and zlib. htslib is vendored — see `lib/` below.
`make check` runs the regression suite.

The Makefile selects a platform block by commenting; the default is
`macos-arm`. Only `lib/macos-arm/libhts.a` is committed, because that is the
only platform it has been built on. For another platform run:

    ./lib/build_htslib.sh linux        # or osx, macos-arm

which fetches htslib, configures it without bzip2/lzma/libcurl/libdeflate
(CRAM features divstats does not use, so nothing beyond zlib is needed) and
writes `lib/<platform>/libhts.a`.

## Input

One chromosome per run, sorted by position, biallelic sites.

- `--vcf FILE` — VCF, bgzipped VCF, or BCF. The format is detected from the
  file; there is no separate flag. `GT` may sit anywhere in `FORMAT`.
- `--tped FILE` — transposed PLINK. Note the missing-allele encodings differ
  between the two readers: VCF uses `.`, TPED uses `-9`.
- `--hemi` — one haplotype per sample instead of two.

Multi-chromosome and unsorted input are refused with a diagnostic rather than
analysed incorrectly.

## Windows

    --sites --winsize 100 --winstep 100     # 100 SNPs, tiling
    --bp --winsize 10000 --winstep 1000     # 10 kb, sliding

### Genetic distance

`--ehh` places its sub-windows by physical position, and does **not** need a
genetic map — earlier versions required `--map` or `--pmap` for any EHH
calculation and then ignored the map entirely. `--ehh-cm` is the placement that
uses it: widths are in the map's own units (normally cM), so they are
comparable between regions of differing recombination rate, and a fixed genetic
width covers fewer base pairs inside a recombination hotspot.

The two centre their sub-windows differently: `--ehh` on the window's
coordinate midpoint, `--ehh-cm` on the genetic midpoint of the window's first
and last SNP. `--ehh-cm` is not computed per partition.

`--partition a b c` splits each window into sub-windows. The values are
percentages under `--sites` and base pairs under `--bp`, and must sum to 100
or to `--winsize` respectively.

## Statistics

| flag | column(s) | meaning |
|---|---|---|
| `--pi` | `pi` | mean pairwise differences per window |
| `--s` | `S` | segregating sites (fractional when subsampling) |
| `--d` | `D` | Tajima's *D*; `nan` where `S` < 1, see below |
| `--h` | `H` | Fay & Wu's *H* |
| `--pik k [k…]` | `pik` | π among the *k* most frequent haplotypes |
| `--ehh w [w…]` | `ehh<w>` | EHH in sub-windows of *w* SNPs (`--sites`) or *w* bp (`--bp`) |
| `--ehhk k w [w…]` | `ehhk<k>_<w>` | EHH among the *k* most frequent haplotypes |
| `--ehh-cm w [w…]` | `ehhcm_<w>` | EHH in sub-windows *w* wide in genetic distance; requires `--map` |

With `--partition`, each statistic also appears per partition, suffixed `_A`,
`_B`, … in partition order.

### Tajima's *D* requires at least one segregating site

`D` is written as the na-string in any window whose `S` is below 1.

Subsampling makes `S` an *expected* count — the number of sites still
polymorphic after projecting to the run's sample size — so a window can report
`0 < S < 1`, which a whole-number count cannot produce. `D` is not merely
noisy there. Its numerator is linear in `S` while its denominator grows as
`sqrt(S)`, so the value is systematically shrunk toward zero by

    D(S) / D(S=1) = sqrt(S) * sqrt(e1 / (e1 - e2*(1 - S)))

exactly — the allele-frequency terms cancel, so the attenuation depends only
on `S` and the sample size. A window holding a third of an expected
segregating site would otherwise report a `D` at about 58% of what that one
site gives, in the same column as windows carrying hundreds of sites.

`S`, `pi` and `H` are still reported for these windows, so they remain
distinguishable from empty ones — filter on `S` if you want them back. Nothing
changes without subsampling, where `S` is a whole number and `S < 1` means
`S = 0`. On projected data the rule is rare: roughly 0.1% of windows at 5%
missing genotypes and 1.1% at 25%, and none once windows average ten or more
segregating sites.

`H` here is the unnormalised *π* − *θ*<sub>H</sub> and `pi` does not involve
`S`, so neither is affected by this rule.

## SweepFinder2 export

    --sweepfinder

Writes `<out>.sweepfinder.out` (`position`, `x`, `n`, `folded`) and exits. `n`
is the number of haplotypes observed at each site, so it varies with local
missingness. `folded` is always 0, meaning the ALT allele is assumed to be
derived — divstats has no outgroup and cannot check this. Sites with no called
genotype are omitted.

## Sample size

Every window's spectrum is projected to **one sample size shared by the whole
run**, reported in the `nhaps` column. Before 2.0.0 each window was projected
to its own minimum observed sample size, so `n` varied with local missingness
and π, S, *D* and *H* were not comparable between windows — and `n` was not in
the output, so nothing downstream could detect it.

The default is the largest `n` reachable by **99.9% of sites**. Taking the
largest reachable by *every* site instead lets the single worst-covered site
set `n` for the whole run — on data with a coverage tail that can cost more
than half the sample for the sake of a handful of sites. Up to 0.1% of sites
are therefore excluded, and the exclude-nothing value is printed so it can be
restored exactly:

    Projecting every window to n = 6 haplotypes (1998 of 2000 sites usable, 99.90%).
      n = 4 would exclude no site; --target-n sets either explicitly.

Files with fewer than 1000 sites are unaffected — 0.1% rounds to zero sites, so
nothing is dropped on a fraction too small to resolve.

`--target-n N` projects to `N` instead and excludes sites observed in fewer
than `N` haplotypes — they cannot be projected upward. Whenever any site is
excluded, including under the default, the number that actually contributed
appears in an `nSNPsUsed` column.

`--window-n-sub` restores the pre-2.0.0 per-window minimum. It maximises each
window's `n` in isolation, but the windows are then not comparable.
`--const-n-sub` is accepted and does nothing; it is the default.

With `--no-sfs-sub` nothing is projected and `nhaps` reports the full sample
size. For statistics that do not use the spectrum (EHH and its variants),
`nhaps` is likewise the full sample size, since those use every haplotype.

## Missing data

Per-site sample size varies with missingness, so windows are not directly
comparable unless the sample size is held fixed. By default divstats projects
each window's spectrum down to that window's minimum sample size
(hypergeometric subsampling).

- default — project to the per-window minimum
- `--const-n-sub` — project to a single global sample size instead, so all
  windows share one *n*
- `--no-sfs-sub` — no projection; each site contributes at its own sample size

Only `pi`, `S`, `D` and `H` are affected; EHH and `pik` are not.

## Output

Tab-separated, to `<out>.divstats.out`:

    chr  start  end  nbps  nSNPs  <statistic columns…>

`start`/`end` are the window's first and last SNP positions, `nbps` its span,
`nSNPs` the sites it contains. Undefined windows are written as `-999`.

## Reproducibility

Statistics changed in 2.0.0 for any input with missing genotypes — see
`CHANGELOG.md`. Do not pool output across that boundary. `divstats --version`
is not yet implemented; the version prints in the startup banner on stderr.

## License

GPLv3. See `LICENSE.md`.
