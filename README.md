# divstats

Window-based diversity statistics for genetic data analysis.

Computes nucleotide diversity, segregating sites, Tajima's *D*, Fay & Wu's
*H*, extended haplotype homozygosity and *k*-haplotype variants of π and EHH,
in windows defined either by a number of SNPs or by physical distance.

## Version

    divstats --version        # prints e.g. 2.0.0 on stdout, exits 0

## Build

    cd src && make

Requires a C++ compiler and zlib. Nothing in the Makefile needs editing — the
platform is detected from `uname`, and the x86 codegen flags are applied only
on x86. Other targets:

    make check      build, then run the regression suite
    make clean      remove objects and the binary
    make install    install to $(PREFIX)/bin   (default /usr/local)
    make info       print the detected platform and flags

Any variable can be overridden on the command line:

    make CXX=clang++ OPT="-O2 -g"
    make install PREFIX=$HOME/.local

### htslib

htslib is vendored as a prebuilt static archive, one per platform, at
`lib/<platform>/libhts.a`. An archive is committed for each platform it has
been built on so far; to see whether yours is covered:

    make info          # prints the detected platform and whether libhts.a is there
    ls lib/*/libhts.a  # what is committed

If it is missing, `make` stops with that path and the command to produce it
rather than a page of undefined symbols:

    ./lib/build_htslib.sh              # platform detected from uname
    ./lib/build_htslib.sh linux 1.21   # or state it explicitly

The script fetches the htslib **release tarball** (the GitHub tag archive
lacks the generated `configure`), builds it without
bzip2/lzma/libcurl/libdeflate — CRAM features divstats does not use, so
nothing beyond zlib is needed — writes `lib/<platform>/libhts.a`, refreshes
`include/htslib/`, and then verifies the archive references none of those
back-ends, failing if it does. Commit the result so nobody on that platform
has to repeat it.

A committed archive is tied to the toolchain and C library it was built
against. If one is present for your platform but linking fails, rebuild it
locally with the same script; the headers are platform-independent and do not
need regenerating.

## Input

One chromosome per run, sorted by position, biallelic sites.

- `--vcf FILE` — VCF, bgzipped VCF, or BCF. The format is detected from the
  file; there is no separate flag. `GT` may sit anywhere in `FORMAT`.
- `--tped FILE` — transposed PLINK. Note the missing-allele encodings differ
  between the two readers: VCF uses `.`, TPED uses `-9`.
- `--hemi` — one haplotype per sample instead of two.

Multi-chromosome and unsorted input are refused with a diagnostic rather than
analysed incorrectly.

## Performance

Parsing dominates a typical run, and it is mostly the cost of turning VCF
text into records — work that does not divide across threads. **BCF input is
substantially faster and produces byte-identical output.** On a 400-haplotype
× 50,000-site file, 500 windows of 100 SNPs, 4 threads:

| input | parse | total |
|---|---|---|
| plain gzip VCF | 0.222 s | 0.348 s |
| bgzipped VCF | 0.222 s | 0.328 s |
| BCF | 0.067 s | 0.203 s |

Convert once with `bcftools view -O b -o data.bcf data.vcf.gz`.

`--pik` cost per window scales with the square of the number of *distinct*
haplotypes in that window, so it is far cheaper on data with missing
genotypes (which exclude a haplotype from the spectrum) than on complete
data, and it grows quickly with window width.

Genotypes are held two bits per site rather than one byte, so the matrix is
`nhaps × nloci / 4` bytes: 4.8 MB for 400 haplotypes × 50,000 sites, and
4.7 GB rather than 18.6 GB for 2,000 haplotypes × 10M sites. `--threads`
covers both the statistics and htslib's decoder, but the decoder only
parallelises BGZF inflation — it does nothing for a plain-gzip VCF, and even
for BGZF it is a few percent, because the text parse is serial.

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

`--partition a b c` splits each window into sub-windows. The values are in
the same units as `--winsize` — SNP counts under `--sites`, base pairs under
`--bp` — and must sum to `--winsize` in both modes. Sub-windows tile the
window with no gap or overlap, in order, and the statistic columns are
suffixed `_A`, `_B`, … in that order. The span each suffix covers is written
to `<out>.divstats.log`.

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
`_B`, … in partition order. What each suffix covers is recorded in
`<out>.divstats.log` — it is not recoverable from the table.

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

## Run log

Every run writes `<out>.divstats.log` alongside the table: version, timestamp,
the full command line, input and window parameters, the subsampling target,
and — when `--partition` is used — a legend giving the span each column suffix
covers:

    partitions: --partition 25 50 25
      Column suffixes _A, _B, ... index the partitions below, in order.
      Spans are relative to the start of each window, not to the chromosome.
      _A  SNPs 1 to 25 of each window  (width 25)
      _B  SNPs 26 to 75 of each window  (width 50)
      _C  SNPs 76 to 100 of each window  (width 25)

This is a separate file rather than `#` comment lines in the table so that
`pd.read_csv(path, sep='\t')` keeps working with no extra arguments.

## Precision

    --precision 12

Significant digits per statistic. The default, 6, is what every previous
version emitted, so output is unchanged unless you ask. Accepts 1–17; a
double carries at most 17.

## Getting help

    divstats --help        # grouped option list with examples; exits 0
    divstats --version     # version to stdout; exits 0
    divstats               # usage; exits 1

Options are grouped under *Input*, *Windows*, *Statistics*, *Sample size and
missing data*, *Output* and *Other*, in that order, with worked examples after
the list.

## Progress

Runs print a line to stderr at each 10% of windows completed, but only once a
run has been going for two seconds — short runs stay quiet. Output goes to
stdout and the output file, so progress can be discarded with `2>/dev/null`
without losing anything.

## Reproducibility

Statistics changed in 2.0.0 for any input with missing genotypes — see
`CHANGELOG.md`. Do not pool output across that boundary. `divstats --version`
is not yet implemented; the version prints in the startup banner on stderr.

## License

GPLv3. See `LICENSE.md`.
