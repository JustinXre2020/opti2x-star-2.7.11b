#!/usr/bin/env bash
# Builds and runs the opti2x differential unit tests against an upstream STAR 2.7.11b source tree.
# Usage: UPSTREAM_SRC=/path/to/STAR-2.7.11b/source opti2x/tests/run_unit_tests.sh [test_name ...]
# Tests that need the rest of STAR link the objects of a built tree (STAR_OBJ_DIR, default the candidate source/).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
SRC=${CAND_SRC:-$HERE/../../source}
: "${UPSTREAM_SRC:?set UPSTREAM_SRC to an upstream STAR 2.7.11b source/ directory}"
CXX=${CXX:-g++-16}
BUILD=${BUILD_DIR:-$HERE/build}
OBJ_DIR=${STAR_OBJ_DIR:-$SRC}
mkdir -p "$BUILD"
CXXFLAGS=(-std=c++11 -O2 -Wall -fopenmp -I"$SRC" -I"$SRC/htslib")

ALL_TESTS="test_extendAlign_diff test_winBinReset test_sjMotif test_stitchAlignToTranscript_diff test_sequenceFuns_diff test_saSearch_diff test_fileMapping test_sjdbInsertSA_diff"
# test name -> extra translation units it needs from $SRC
units() {
  case "$1" in
    test_extendAlign_diff) echo "extendAlign.cpp Transcript.cpp" ;;
    test_winBinReset) echo "" ;;
    test_sjMotif) echo "sjMotif.cpp" ;;
    test_stitchAlignToTranscript_diff) echo "stitchAlignToTranscript.cpp sjMotif.cpp" ;;
    test_sequenceFuns_diff) echo "SequenceFuns.cpp" ;;
    test_saSearch_diff) echo "saSearchPiece.cpp" ;;
    test_fileMapping) echo "fileMapping.cpp" ;;
    test_sjdbInsertSA_diff) echo "sjdbInsertSA.cpp PackedArray.cpp" ;;
    *) echo "unknown test $1" >&2; exit 2 ;;
  esac
}
# tests that link the other STAR objects (the units above are compiled from $SRC instead)
needs_star() { [ "$1" = test_stitchAlignToTranscript_diff ] || [ "$1" = test_saSearch_diff ]; }
# the Makefile OBJECTS list (a few entries are .cpp files compiled at link time), minus STAR.o and the test's own units
star_objs() {
  local skip=" STAR.o $(units "$1" | sed 's/\.cpp/.o/g') " f b
  for f in $(sed -n '/^OBJECTS *=/,/[^\\]$/p' "$SRC/Makefile" | sed 's/^OBJECTS *=//; s/\\$//'); do
    b=${f%.o}; b=${b%.cpp}; b=${b%.c}
    case "$skip" in *" $b.o "*) continue ;; esac
    case "$f" in *.o) echo "$OBJ_DIR/$f" ;; *) echo "$SRC/$f" ;; esac
  done
}
DEFS=(-DUPSTREAM_EXTENDALIGN="\"$UPSTREAM_SRC/extendAlign.cpp\"" -DUPSTREAM_STITCH="\"$UPSTREAM_SRC/stitchAlignToTranscript.cpp\""
      -DUPSTREAM_SEQUENCEFUNS="\"$UPSTREAM_SRC/SequenceFuns.cpp\"" -DSTAR_BIN="\"${STAR_BIN:-$OBJ_DIR/STAR}\"")

if [ $# -gt 0 ]; then TESTS=("$@"); else TESTS=($ALL_TESTS); fi

rc=0
for t in "${TESTS[@]}"; do
  srcs=()
  for u in $(units "$t"); do srcs+=("$SRC/$u"); done
  link=()
  if needs_star "$t"; then link=($(star_objs "$t") -fopenmp -pthread -lz "$SRC/htslib/libhts.a"); fi
  "$CXX" "${CXXFLAGS[@]}" "${DEFS[@]}" -o "$BUILD/$t" "$HERE/$t.cpp" ${srcs[@]+"${srcs[@]}"} ${link[@]+"${link[@]}"}
  if "$BUILD/$t"; then echo "PASS $t"; else echo "FAIL $t"; rc=1; fi
done
exit $rc
