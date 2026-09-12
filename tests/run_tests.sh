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
SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/divstats-tests.XXXXXX")"
cleanup() { [ "$KEEP" = 1 ] && echo "scratch kept: $SCRATCH" || rm -rf "$SCRATCH"; }
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
        if (a==e) continue
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

define_case sites-basic       table -- --vcf "$C" --sites --winsize 100 --winstep 100 --pi --s --d --h
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
define_case threads-4         table -- --vcf "$C" --sites --winsize 100 --winstep 100 --pi --s --d --h --threads 4
# ===========================================================================

golden_for() {   # thread-invariance shares sites-basic's golden
  case "$1" in threads-4) echo "sites-basic" ;; *) echo "$1" ;; esac
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
    [ "$name" = "threads-4" ] && { printf "  ----  %-18s shares sites-basic golden\n" "$name"; return; }
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
