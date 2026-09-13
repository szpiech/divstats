#!/usr/bin/env bash
# build_htslib.sh -- produce the vendored libhts.a for one platform
#
#   ./lib/build_htslib.sh [platform] [version]
#
# platform defaults to a guess from uname (macos-arm, osx, linux); version
# defaults to the value below. Writes lib/<platform>/libhts.a and refreshes
# include/htslib/*.h.
#
# Only lib/macos-arm/libhts.a is committed, because that is the only platform
# the library has been built on here. Run this on a Linux box to populate
# lib/linux, and so on. The headers are platform-independent, so they only
# need generating once.
#
# The configuration is deliberately minimal: bzip2, lzma, libcurl and
# libdeflate are CRAM features that divstats does not use, and disabling them
# means libhts.a needs nothing beyond the -lz that divstats already links.
# Check that after building with:
#
#   nm -u lib/<platform>/libhts.a | grep -E '_(BZ2|lzma|curl|libdeflate)'
#
# which should print nothing.
set -euo pipefail

VERSION="${2:-1.21}"

if [ -n "${1:-}" ]; then
  PLATFORM="$1"
else
  case "$(uname -s)-$(uname -m)" in
    Darwin-arm64)  PLATFORM=macos-arm ;;
    Darwin-x86_64) PLATFORM=osx ;;
    Linux-*)       PLATFORM=linux ;;
    *) echo "error: cannot guess platform from $(uname -s)-$(uname -m); pass it explicitly" >&2
       exit 2 ;;
  esac
fi

HERE="$(cd "$(dirname "$0")" && pwd -P)"
ROOT="$(cd "$HERE/.." && pwd -P)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/htslib-build.XXXXXX")"
trap 'rm -rf -- "$WORK"' EXIT

echo "building htslib $VERSION for $PLATFORM in $WORK"
cd "$WORK"

# Use the RELEASE tarball, not the git tag archive: the tag archive does not
# contain the generated ./configure, so it needs autoconf installed.
URL="https://github.com/samtools/htslib/releases/download/${VERSION}/htslib-${VERSION}.tar.bz2"
curl -fsSL -o hts.tar.bz2 "$URL"
tar xjf hts.tar.bz2
cd "htslib-${VERSION}"

./configure --disable-libcurl --disable-bz2 --disable-lzma \
            --without-libdeflate --disable-plugins --disable-gcs --disable-s3
make -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)" lib-static

mkdir -p "$ROOT/lib/$PLATFORM" "$ROOT/include/htslib"
cp -p libhts.a "$ROOT/lib/$PLATFORM/libhts.a"
cp -p htslib/*.h "$ROOT/include/htslib/"

echo "wrote lib/$PLATFORM/libhts.a ($(wc -c < libhts.a) bytes) and include/htslib/"
echo "checking for unwanted external dependencies:"
if nm -u libhts.a 2>/dev/null | grep -qE '_(BZ2|lzma|curl|libdeflate)'; then
  echo "  WARNING: libhts.a references bz2/lzma/curl/libdeflate symbols;"
  echo "           divstats links only -lz and will fail at link time."
  exit 1
fi
echo "  none -- only zlib and libc are needed."
