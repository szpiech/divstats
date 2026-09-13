#!/usr/bin/env bash
# run_tests.sh -- golden-file regression suite for divstats
#
#   ./tests/run_tests.sh              run every case
#   ./tests/run_tests.sh sites-basic  run named case(s)
#   ./tests/run_tests.sh --list       list case names
#   ./tests/run_tests.sh --regen      overwrite the goldens from the current binary
#   ./tests/run_tests.sh --bin PATH   test a specific binary (default src/divstats)
#   ./tests/run_tests.sh --keep       keep the scratch directory and print its path
#
# Numeric tables are compared COLUMN BY COLUMN: the header must match exactly,
# integer-looking columns must match exactly, floats must agree to DEFAULT_RTOL.
# A hash would tell you that something moved but never which statistic.
#
# Fixtures are derived at run time from tests/data/core.vcf.gz -- a 400-site
# excerpt of test/test.small.vcf, the only input committed to the repo. Cases
# needing the large unversioned files under test/ are skipped when absent.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
BIN="$ROOT/src/divstats"
EXPECTED="$HERE/expected"
CORE="$HERE/data/core.vcf.gz"
DEFAULT_RTOL="1e-6"
REGEN=0
KEEP=0
LIST=0
SELECT=()

while [ $# -gt 0 ]; do
  case "$1" in
    --regen) REGEN=1 ;;
    --keep)  KEEP=1 ;;
    --list)  LIST=1 ;;
    --bin)   shift; BIN="$1" ;;
    --rtol)  shift; DEFAULT_RTOL="$1" ;;
    -h|--help) sed -n '2,20p' "$BASH_SOURCE"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 2 ;;
    *)  SELECT+=("$1") ;;
  esac
  shift
done

# ---------------------------------------------------------------- portability
hash_file() { if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1
              else md5 -q "$1"; fi; }
GZCAT="gzip -dc"

PASS=0; FAIL=0; SKIP=0
FAILED_CASES=()

# ------------------------------------------------------------------ scratch
# mktemp can fail -- e.g. TMPDIR set to a relative or nonexistent path -- and a
# bare "rm -rf $SCRATCH" in the exit trap then targets whatever the variable
# expanded to. Validate before anything is registered for deletion, and have
# cleanup refuse any path that is not a scratch directory of ours.
SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/divstats-tests.XXXXXX" 2>/dev/null)" || SCRATCH=""
if [ -z "$SCRATCH" ] || [ ! -d "$SCRATCH" ]; then
  echo "error: could not create a scratch directory (TMPDIR='${TMPDIR:-unset}')" >&2
  exit 2
fi
SCRATCH="$(cd "$SCRATCH" && pwd -P)"
case "$SCRATCH" in
  */divstats-tests.??????) : ;;
  *) echo "error: unexpected scratch path '$SCRATCH'; refusing to run" >&2; exit 2 ;;
esac

cleanup() {
  case "${SCRATCH:-}" in
    */divstats-tests.??????) : ;;
    *) return 0 ;;                      #never remove anything else
  esac
  [ -d "$SCRATCH" ] || return 0
  if [ "$KEEP" = 1 ]; then
    echo "scratch kept: $SCRATCH"
  else
    rm -rf -- "$SCRATCH"
  fi
}
trap cleanup EXIT

# ------------------------------------------------------------ fixture build
# All transforms are deterministic -- no RNG -- so the goldens are stable.
build_fixtures() {
  $GZCAT "$CORE" > "$SCRATCH/core.vcf"

  # missing-data fixture: on every 3rd variant, blank every 7th sample.
  # Exercises the SFS subsampling path, which no shipped input covers with a
  # small enough file to keep the suite fast.
  awk 'BEGIN{OFS="\t"}
       /^#/ {print; next}
       {v++
        if (v%3==0) for (i=10;i<=NF;i++) if ((i-9)%7==0) $i="./."
        print}' "$SCRATCH/core.vcf" > "$SCRATCH/missing.vcf"

  # hemizygous fixture: keep the first allele of each genotype (cf. --hemi)
  awk 'BEGIN{OFS="\t"}
       /^#/ {print; next}
       {for (i=10;i<=NF;i++) $i=substr($i,1,1); print}' "$SCRATCH/core.vcf" \
       > "$SCRATCH/hemi.vcf"

  # 4-column map: <chr> <locusID> <genetic pos> <physical pos>.
  # Genetic position is a nominal 1 cM/Mb; divstats never reads it (the genetic
  # map is parsed and ignored), but --ehh demands the file unless --pmap is set.
  awk '!/^#/ {printf "%s\tsnp%d\t%.6f\t%d\n", $1, ++n, $2/1000000.0, $2}' \
      "$SCRATCH/core.vcf" > "$SCRATCH/core.map"

  # TPED + map, to cover the --tped reader at all.
  # Note the encoding difference: the VCF reader takes "." for a missing
  # allele, the TPED reader requires "-9" (TPED_MISSING). Neither is
  # documented in --help.
  awk '!/^#/ {printf "%s\tsnp%d\t%.6f\t%d", $1, ++n, $2/1000000.0, $2
              for (i=10;i<=NF;i++) {
                 split($i, g, /[|\/]/)
                 a = (g[1]=="." ? "-9" : g[1]); b = (g[2]=="." ? "-9" : g[2])
                 printf "\t%s\t%s", a, b }
              printf "\n"}' "$SCRATCH/core.vcf" > "$SCRATCH/core.tped"

  gzip -c "$SCRATCH/missing.vcf" > "$SCRATCH/missing.vcf.gz"
  gzip -c "$SCRATCH/hemi.vcf"    > "$SCRATCH/hemi.vcf.gz"

  # Malformed inputs, for the reject cases. Both used to be accepted silently
  # and produce wrong windows.
  #   multichr : the second half of the records relabelled to another contig
  #              already declared in the header
  #   unsorted : the same records in reverse position order
  awk 'BEGIN{OFS="\t"} /^#/{print; next}
       {n++; if (n > 200) $1="chr10"; print}' "$SCRATCH/core.vcf" \
     | gzip -c > "$SCRATCH/multichr.vcf.gz"
  awk '/^#/{print; next} {rec[++n]=$0}
       END{for(i=n;i>=1;i--) print rec[i]}' "$SCRATCH/core.vcf" \
     | gzip -c > "$SCRATCH/unsorted.vcf.gz"
}

# --------------------------------------------------------------- comparators
# compare_table <expected.tsv.gz> <actual.tsv> <rtol> [tol_col_regex] [max_frac]
# Splits on any whitespace so it works both before and after the header
# delimiter fix (the header used to be space-joined while rows were tabs).
compare_table() {
  local exp="$1" act="$2" rtol="${3:-$DEFAULT_RTOL}" tolcol="${4:-}" maxfrac="${5:-0}"
  $GZCAT "$exp" > "$SCRATCH/.exp" 2>/dev/null || { echo "    no golden: $exp"; return 1; }
  awk -v rtol="$rtol" -v tolcol="$tolcol" -v maxfrac="$maxfrac" '
    function abs(x){ return x<0 ? -x : x }
    function isnum(s){ return s ~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/ }
    function isint(s){ return s ~ /^[+-]?[0-9]+$/ }
    # nan/inf must be handled EXPLICITLY, before any use of == or <, and the
    # comparison must be forced to strings. awk parses "nan" with strtod and
    # then compares numerically through a three-way (a<b ? -1 : a>b ? 1 : 0);
    # both tests are false for NaN, so "nan" compares EQUAL to every number.
    # Without this the suite silently passes a build that emits nan for every
    # statistic -- which is exactly what a cohort past the nCk overflow
    # threshold produces, and exactly what the log-space weights fix.
    function nonfinite(s){ return tolower(s) ~ /^[+-]?(nan|inf|infinity)$/ }
    NR==FNR {
      n=split($0, f, /[ \t]+/); erows++
      for (i=1;i<=n;i++) E[erows,i]=f[i]
      ecols[erows]=n
      next
    }
    {
      n=split($0, f, /[ \t]+/); arows++
      if (arows==1) {                                  # header: exact, in order
        if (n!=ecols[1]) { printf "    header: %d fields, golden has %d\n", n, ecols[1]; bad=1 }
        for (i=1;i<=n && i<=ecols[1];i++) {
          name[i]=f[i]
          if (f[i]!=E[1,i]) { printf "    header field %d: %s, golden %s\n", i, f[i], E[1,i]; bad=1 }
        }
        next
      }
      if (arows>erows) { extra++; next }
      for (i=1;i<=n;i++) {
        a=f[i]; e=E[arows,i]
        if (nonfinite(a) || nonfinite(e)) {
          # identical tokens (nan vs nan) match; anything else is a difference
          if (("" a) != ("" e)) { diff[i]++; md[i]=-1 }
          continue
        }
        if (("" a) == ("" e)) continue
        if (isnum(a) && isnum(e)) {
          if (isint(a) && isint(e)) { diff[i]++; if (1>md[i]) md[i]=1; continue }
          d = abs(a-e) / (abs(e)>1e-300 ? abs(e) : 1)
          if (d > rtol) { diff[i]++; if (d>md[i]) md[i]=d }
        } else { diff[i]++; md[i]=-1 }
      }
    }
    END {
      if (arows!=erows) { printf "    row count: %d, golden has %d\n", arows-0, erows-0; bad=1 }
      ndata = (erows>1 ? erows-1 : 1)
      for (i in diff) {
        frac = diff[i]/ndata
        tolerated = (tolcol!="" && name[i] ~ tolcol && frac <= maxfrac+1e-12)
        if (md[i] < 0)
          printf "    column %-14s %4d/%d rows differ, non-numeric (nan/inf or text)%s\n",
                 name[i], diff[i], ndata, (tolerated ? "  [tolerated]" : "")
        else
          printf "    column %-14s %4d/%d rows differ, max rel dev %.3g%s\n",
                 name[i], diff[i], ndata, md[i], (tolerated ? "  [tolerated]" : "")
        if (!tolerated) bad=1
      }
      exit (bad ? 1 : 0)
    }
  ' "$SCRATCH/.exp" "$act"
}

compare_hash() {
  local exp="$1" act="$2"
  local want have
  want="$(cat "$exp" 2>/dev/null)" || { echo "    no golden: $exp"; return 1; }
  have="$(hash_file "$act")"
  [ "$want" = "$have" ] && return 0
  echo "    hash $have, golden $want"
  return 1
}

# ------------------------------------------------------------------ case run
# run_case <name> <output-kind> <golden-basename> -- <divstats args...>
# output-kind: table | stdout | smoke:<expected-row-count>
declare -a CASE_NAMES=()
declare -a CASE_KIND=()
declare -a CASE_ARGS=()
declare -a CASE_TOLCOL=()
declare -a CASE_MAXFRAC=()

define_case() {           # define_case <name> <kind> [tolcol] [maxfrac] -- args...
  local name="$1" kind="$2"; shift 2
  local tolcol="" maxfrac=0
  if [ "${1:-}" != "--" ]; then tolcol="$1"; maxfrac="$2"; shift 2; fi
  shift                   # the --
  CASE_NAMES+=("$name"); CASE_KIND+=("$kind")
  CASE_TOLCOL+=("$tolcol"); CASE_MAXFRAC+=("$maxfrac")
  CASE_ARGS+=("$*")
}

# ============================== CASE LIST ==================================
# Covers: both window units (--sites/--bp), tiling and sliding windows, every
# statistic, all three missing-data modes, partitions in both window units,
# both EHH map sources, both input formats, --hemi, the SweepFinder export,
# and thread invariance.
C="$SCRATCH/core.vcf"
M="$SCRATCH/missing.vcf.gz"
H="$SCRATCH/hemi.vcf.gz"

# Hand-computable fixture: 4 diploid samples, 6 sites, one missing genotype at
# site 4, so the window's target sample size is H = 6 while five sites carry
# n = 8. The hypergeometric projection is small enough to evaluate exactly:
#
#   projected SFS  0.5, 1.714286, 1.071429, 2.071429, 0.535714, 0.107143, 0
#   pi = 2.707143   S = 5.5   Tajima's D = 1.322279 (at S = 5.5: 0.701735)
#
# Before the population-size fix divstats reported pi = 2.6, matching a
# projection from n+1 chromosomes. This case is the arithmetic check on the
# subsampler; the numbers above are independent of the implementation.
define_case subsample-toy     table -- --vcf "$HERE/data/subsample-toy.vcf.gz" --sites --winsize 6 --winstep 6 --pi --s --d
# 1200 haplotypes, 5% missing. The old weights formed each binomial
# coefficient separately, and it is the PRODUCT nCk(i,j)*nCk(n-i,H-j) that
# overflows first -- two mid-range coefficients near C(600,300) ~ 1e179
# multiply to inf long before C(n,H) itself is large. inf/inf is nan, so every
# SFS statistic came back nan on this fixture.
#
# Measured with the pre-fix build (pi for a single 20-site window):
#     1040 haplotypes,  1% missing -> nan
#     1040 haplotypes,  5% missing -> 4.96285
#     1040 haplotypes, 10% missing -> 4.94862
#     1200 haplotypes,  1% missing -> nan
#     1200 haplotypes,  5% missing -> nan
# So the trigger is NOT monotone in the missing-data rate: what matters is
# where H (the per-window minimum sample size) lands relative to the per-site
# sample sizes, which decides how close j and H-j get to the middle of their
# ranges. This fixture was chosen empirically, not derived.
#
# Threshold is platform-dependent -- with 80-bit long double it is far higher
# -- so on x86-64 this case guards the arithmetic without necessarily
# reproducing the original failure.
define_case large-n-subsample table -- --vcf "$HERE/data/large-n.vcf.gz" --sites --winsize 20 --winstep 20 --pi --s --d --h
# Tajima's D with no segregating sites in the window: 119 of the 200 windows
# here have S = 0, where the numerator and the variance estimate are both 0.
# This pins the result at the na-string rather than 0, inf, or whatever
# 0.0/0.0 happens to produce. Also exercises --na-string itself.
define_case tinywin-undef    table -- --vcf "$C" --sites --winsize 2 --winstep 2 --pi --s --d --h
define_case na-string-NA     table -- --vcf "$C" --sites --winsize 10 --winstep 10 --pik 4 --na-string NA
define_case sites-basic       table -- --vcf "$C" --sites --winsize 100 --winstep 100 --pi --s --d --h
# Same arguments as sites-basic against the same variants stored as BCF, and
# it shares sites-basic's golden -- so it asserts that the binary format gives
# identical statistics to the text one, not merely that BCF parses.
define_case bcf-basic         table -- --vcf "$HERE/data/core.bcf" --sites --winsize 100 --winstep 100 --pi --s --d --h
define_case sites-pi-only     table -- --vcf "$C" --sites --winsize 100 --winstep 100 --pi
define_case sites-sliding     table -- --vcf "$C" --sites --winsize 100 --winstep 25  --pi --s --d --h
define_case bp-basic          table -- --vcf "$C" --bp --winsize 200000 --winstep 200000 --pi --s --d --h
define_case bp-sliding        table -- --vcf "$C" --bp --winsize 200000 --winstep 50000 --pi --s
define_case missing-default   table -- --vcf "$M" --sites --winsize 100 --winstep 100 --pi --s --d --h
define_case missing-nosub     table -- --vcf "$M" --sites --winsize 100 --winstep 100 --pi --s --d --h --no-sfs-sub
define_case missing-constn    table -- --vcf "$M" --sites --winsize 100 --winstep 100 --pi --s --d --h --const-n-sub
define_case ehh-pmap          table -- --vcf "$C" --sites --winsize 100 --winstep 100 --ehh 20 50 --pmap
define_case ehh-mapfile       table -- --vcf "$C" --sites --winsize 100 --winstep 100 --ehh 20 50 --map "$SCRATCH/core.map"
define_case ehhk              table -- --vcf "$C" --sites --winsize 100 --winstep 100 --ehh 20 50 --ehhk 2 4 --pmap
define_case pik               table -- --vcf "$C" --sites --winsize 50  --winstep 50  --pik 2
# These two reach the pi_k2 path that used to read past the end of
# hfs->sortedCount (review finding B2). The over-read needs
# distinct-counts < k <= distinct-haplotypes, which is the small-window /
# moderate-k regime -- a large k on a wide window breaks out on the first
# iteration and never overflows. Confirmed with ASan on the pre-fix build.
define_case pik-narrow-window table -- --vcf "$C" --sites --winsize 10  --winstep 10  --pik 4
define_case pik-multi-k       table -- --vcf "$C" --sites --winsize 50  --winstep 50  --pik 2 3 4 5 6 7 8
# The three cases above only ever reach the tie-averaging branch with t=2,
# m=1 -- the trivial tie, where one haplotype is drawn from a class of two.
# Confirmed by instrumenting the branch with a counter. These two reach
# t=4, m=2 and t=4, m=3, which is where the closed form's third term
# (tied x tied pairs, coefficient m(m-1)/(t(t-1))) is actually nonzero.
# Without them the closed form's interesting case would be uncovered.
define_case pik-tie-m2        table -- --vcf "$C" --sites --winsize 200 --winstep 200 --pik 2
define_case pik-tie-m3        table -- --vcf "$C" --sites --winsize 200 --winstep 200 --pik 3
define_case partition-sites   table -- --vcf "$C" --sites --winsize 100 --winstep 100 --partition 25 50 25 --pi --s --d
# --ehh-part was unreachable before the validation order was fixed (B7):
# it was tested against DO_PARTITION before --partition had been parsed.
define_case ehh-part          table -- --vcf "$C" --sites --winsize 100 --winstep 100 --partition 25 50 25 --ehh 20 --ehh-part --pi --pmap
define_case partition-bp      table -- --vcf "$C" --bp --winsize 200000 --winstep 200000 --partition 50000 100000 50000 --pi --s
define_case tped-basic        table -- --tped "$SCRATCH/core.tped" --sites --winsize 100 --winstep 100 --pi --s --d --h
define_case hemi              table -- --vcf "$H" --hemi --sites --winsize 100 --winstep 100 --pi --s --d --h
define_case sweepfinder       stdout -- --vcf "$C" --sites --winsize 100 --winstep 100 --sweepfinder
# Thread invariance: compared against the SINGLE-THREADED sites-basic golden,
# so this asserts a property (results independent of --threads), not just
# self-consistency.
# Reject cases: the argument after the kind is a regex the diagnostic must
# match. Both of these inputs used to be accepted silently -- the
# multi-chromosome file was merged into one coordinate space and labelled with
# the last chromosome seen, and the unsorted file produced windows containing
# a subset of the sites with -999 statistics in the first one.
define_case reject-multichr   reject "more than one chromosome" 0 -- --vcf "$SCRATCH/multichr.vcf.gz" --bp --winsize 5000 --winstep 5000 --pi
define_case reject-unsorted   reject "not sorted by position"   0 -- --vcf "$SCRATCH/unsorted.vcf.gz" --bp --winsize 5000 --winstep 5000 --pi
define_case threads-4         table -- --vcf "$C" --sites --winsize 100 --winstep 100 --pi --s --d --h --threads 4
# ===========================================================================

golden_for() {   # thread-invariance shares sites-basic's golden
  case "$1" in threads-4|bcf-basic) echo "sites-basic" ;; *) echo "$1" ;; esac
}

run_one() {
  local idx="$1" name="${CASE_NAMES[$1]}" kind="${CASE_KIND[$1]}"
  local args="${CASE_ARGS[$1]}" tolcol="${CASE_TOLCOL[$1]}" maxfrac="${CASE_MAXFRAC[$1]}"
  local gold; gold="$(golden_for "$name")"
  local out="$SCRATCH/$name"
  local rc

  case "$kind" in
    stdout)
      # shellcheck disable=SC2086
      $BIN $args --out "$out" > "$out.stdout" 2>"$out.err"; rc=$?
      ;;
    *)
      # shellcheck disable=SC2086
      $BIN $args --out "$out" > /dev/null 2>"$out.err"; rc=$?
      ;;
  esac

  # A reject case asserts that bad input is REFUSED, which is the whole point
  # of the multi-chromosome and sort-order checks -- there is no output table
  # to compare, so the assertion is the exit status plus a diagnostic that
  # names the problem. Without this kind, a regression that silently went back
  # to accepting these files would look like a pass.
  if [ "$kind" = "reject" ]; then
    if [ $rc -eq 0 ]; then
      printf "  FAIL  %-18s expected rejection, exited 0\n" "$name"
      FAIL=$((FAIL+1)); FAILED_CASES+=("$name"); return
    fi
    if [ $rc -ge 128 ]; then
      printf "  FAIL  %-18s died on signal (exit %d), expected a clean refusal\n" "$name" $rc
      sed 's/^/    /' "$out.err" | tail -2
      FAIL=$((FAIL+1)); FAILED_CASES+=("$name"); return
    fi
    if ! grep -qE "$tolcol" "$out.err"; then
      printf "  FAIL  %-18s refused, but stderr does not match /%s/\n" "$name" "$tolcol"
      sed 's/^/    /' "$out.err" | head -2
      FAIL=$((FAIL+1)); FAILED_CASES+=("$name"); return
    fi
    printf "  ok    %-18s refused cleanly (exit %d)\n" "$name" $rc
    PASS=$((PASS+1)); return
  fi

  if [ $rc -ne 0 ]; then
    printf "  FAIL  %-18s exited %d\n" "$name" $rc
    sed 's/^/    /' "$out.err" | tail -3
    FAIL=$((FAIL+1)); FAILED_CASES+=("$name"); return
  fi

  local produced golden_path
  case "$kind" in
    stdout) produced="$out.stdout";           golden_path="$EXPECTED/$gold.md5" ;;
    *)      produced="$out.divstats.out";     golden_path="$EXPECTED/$gold.tsv.gz" ;;
  esac

  if [ ! -s "$produced" ]; then
    printf "  FAIL  %-18s produced no output\n" "$name"
    FAIL=$((FAIL+1)); FAILED_CASES+=("$name"); return
  fi

  if [ "$REGEN" = 1 ]; then
    case "$name" in
      threads-4|bcf-basic) printf "  ----  %-18s shares sites-basic golden\n" "$name"; return ;;
    esac
    mkdir -p "$EXPECTED"
    case "$kind" in
      stdout) hash_file "$produced" > "$golden_path" ;;
      *)      gzip -9c "$produced" > "$golden_path" ;;
    esac
    printf "  REGEN %-18s %s\n" "$name" "$(basename "$golden_path")"
    return
  fi

  local msg
  case "$kind" in
    stdout) msg="$(compare_hash  "$golden_path" "$produced")" ; rc=$? ;;
    *)      msg="$(compare_table "$golden_path" "$produced" "$DEFAULT_RTOL" "$tolcol" "$maxfrac")" ; rc=$? ;;
  esac
  if [ $rc -eq 0 ]; then
    printf "  ok    %-18s\n" "$name"; PASS=$((PASS+1))
  else
    printf "  FAIL  %-18s\n" "$name"; [ -n "$msg" ] && echo "$msg"
    FAIL=$((FAIL+1)); FAILED_CASES+=("$name")
  fi
}

# ---------------------------------------------------------------------- main
if [ "$LIST" = 1 ]; then printf '%s\n' "${CASE_NAMES[@]}"; exit 0; fi

if [ ! -x "$BIN" ]; then echo "no binary at $BIN (run make first)" >&2; exit 2; fi
if [ ! -f "$CORE" ]; then echo "missing fixture $CORE" >&2; exit 2; fi

echo "divstats regression suite"
echo "  binary : $BIN"
echo "  rtol   : $DEFAULT_RTOL"
build_fixtures

for i in "${!CASE_NAMES[@]}"; do
  if [ ${#SELECT[@]} -gt 0 ]; then
    keep=0; for s in "${SELECT[@]}"; do [ "$s" = "${CASE_NAMES[$i]}" ] && keep=1; done
    [ $keep = 1 ] || continue
  fi
  run_one "$i"
done

echo "  ---"
printf "  %d passed, %d failed, %d skipped\n" "$PASS" "$FAIL" "$SKIP"
if [ "$FAIL" -gt 0 ]; then printf "  failed: %s\n" "${FAILED_CASES[*]}"; exit 1; fi
exit 0
