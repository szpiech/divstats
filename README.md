# divstats

Window-based diversity statistics for genetic data analysis.

Computes nucleotide diversity, segregating sites, Tajima's *D*, Fay & Wu's
*H*, extended haplotype homozygosity and *k*-haplotype variants of π and EHH,
in windows defined either by a number of SNPs or by physical distance.

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

`--partition a b c` splits each window into sub-windows by percentage.

## Statistics

| flag | column(s) | meaning |
|---|---|---|
| `--pi` | `pi` | mean pairwise differences per window |
| `--s` | `S` | segregating sites (fractional when subsampling) |
| `--d` | `D` | Tajima's *D* |
| `--h` | `H` | Fay & Wu's *H* |
| `--pik k [k…]` | `pik` | π among the *k* most frequent haplotypes |
| `--ehh w [w…]` | `ehh<w>` | EHH in sub-windows of *w* SNPs |
| `--ehhk k w [w…]` | `ehhk<k>_<w>` | EHH among the *k* most frequent haplotypes |

With `--partition`, each statistic also appears per partition, suffixed `_A`,
`_B`, … in partition order.

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
