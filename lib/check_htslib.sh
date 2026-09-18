#!/usr/bin/env bash
# check_htslib.sh -- verify a vendored libhts.a is usable by divstats
#
#   ./lib/check_htslib.sh lib/linux/libhts.a
#
# Two properties, both of which have actually been violated:
#
#   1. Nothing needed beyond zlib and libc. build_htslib.sh configures
#      bzip2, lzma, libcurl and libdeflate away; if that regresses, divstats
#      fails at link time with undefined symbols, because it links only -lz.
#
#   2. Position-independent code. Current Linux toolchains link executables
#      as PIE by default, and a PIE link rejects the absolute 32-bit
#      relocations that non-PIC code emits. The first committed lib/linux
#      archive was built without -fPIC and could not be linked by any default
#      Ubuntu toolchain:
#
#        ld: libhts.a(hts.o): relocation R_X86_64_32 against `.rodata.str1.1'
#            can not be used when making a PIE object; recompile with -fPIE
#
#      That archive was fine for whoever built it and broken for everyone
#      else, which is the failure this script exists to catch.
#
# Run by build_htslib.sh on the archive it has just produced, and by CI on the
# archive that is COMMITTED. Those are different claims, and it was the second
# one that failed.
set -euo pipefail

A="${1:-}"
if [ -z "$A" ]; then
  echo "usage: $(basename "$0") <path to libhts.a>" >&2
  exit 2
fi
if [ ! -f "$A" ]; then
  echo "check_htslib: $A does not exist" >&2
  exit 2
fi

echo "checking $A ($(wc -c < "$A" | tr -d ' ') bytes)"
rc=0

# --- 1. external dependencies -------------------------------------------
printf '  external dependencies: '
if nm -u "$A" 2>/dev/null | grep -qE '_?(BZ2_|lzma_|curl_|libdeflate_)'; then
  echo "FAIL"
  echo "    references bz2/lzma/curl/libdeflate symbols; divstats links only"
  echo "    -lz and will fail at link time. Rebuild with build_htslib.sh,"
  echo "    which configures those back-ends away."
  rc=1
else
  echo "ok -- only zlib and libc"
fi

# --- 2. position-independent code ---------------------------------------
# Any readelf will do. llvm-readelf is accepted deliberately: it reads ELF on
# a non-ELF host, so a macOS user can check a Linux archive -- the case that
# would have caught this before it was committed.
READELF=""
for c in readelf llvm-readelf eu-readelf; do
  if command -v "$c" >/dev/null 2>&1; then READELF="$c"; break; fi
done

printf '  position-independent: '
if [ -z "$READELF" ]; then
  echo "skipped -- no readelf on PATH"
else
  relocs="$("$READELF" -r "$A" 2>/dev/null || true)"
  # Here-strings, not pipes. `printf ... | grep -q` looks equivalent but is
  # not: grep -q exits at the first match and closes the pipe, printf dies of
  # SIGPIPE, and `set -o pipefail` reports the pipeline as failed -- so the
  # test below took the "skipped" branch on an archive that plainly had
  # relocation sections, silently passing the very archive it exists to reject.
  if ! grep -q '^Relocation section' <<< "$relocs"; then
    echo "skipped -- not ELF, or no relocation sections"
  else
    # Relocations in .debug_* sections are absolute even in PIC builds and no
    # linker objects to them, so the section name is tracked rather than
    # grepping the whole listing. On the broken lib/linux archive that
    # distinction was the difference between 3367 real offenders and 90589
    # counting debug info.
    bad="$(awk '
      /^Relocation section/ { sec = $3; gsub(/[^.A-Za-z0-9_]/, "", sec); next }
      sec ~ /debug/ { next }
      $3 == "R_X86_64_32" || $3 == "R_X86_64_32S" || $3 == "R_AARCH64_ABS32" { n++ }
      END { print n + 0 }' <<< "$relocs")"
    if [ "$bad" -gt 0 ]; then
      echo "FAIL"
      echo "    $bad absolute 32-bit relocation(s) outside debug sections, so"
      echo "    this archive cannot be linked into a PIE executable -- the"
      echo "    default on current Linux toolchains. It was built without"
      echo "    -fPIC; rebuild it with build_htslib.sh."
      rc=1
    else
      echo "ok -- no absolute relocations in allocated sections"
    fi
  fi
fi

exit $rc