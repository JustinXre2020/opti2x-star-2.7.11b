#!/usr/bin/env bash
# Junction insertion on the fly must leave every output byte-identical to upstream. With an index built from the same GTF
# (encode config) no suffix array entry changes, so the entries are moved as a block; the second pass of twopassMode
# inserts the junctions the first pass found, which are new, so the entries are rebuilt one by one. With more than one
# thread, the threads share both.
# Usage: STAR_IDX=... STAR_R1=... STAR_R2=... STAR_GTF=... test_sjdbInsert.sh <candidate STAR> <upstream STAR> <workdir> [threads]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
CAND=$1; UPS=$2; WORK=$3; THREADS=${4:-1}
BLOCK="SA entries unchanged by the insertion, moved as a block"

rc=0
"$HERE/oracle_test.sh" "$CAND" "$UPS" "$WORK" "$THREADS" encode twopass || rc=1

log() { echo "$WORK/candidate-$1-t$THREADS/Log.out"; }
if grep -qF -- "$BLOCK" "$(log encode)" 2>/dev/null; then echo "PASS encode log: $BLOCK"; else echo "FAIL encode: $(log encode) lacks: $BLOCK"; rc=1; fi
if grep -q "number of new junctions=[1-9]" "$(log twopass)" 2>/dev/null && ! grep -qF -- "$BLOCK" "$(log twopass)"; then
  echo "PASS twopass log: new junctions, entries rebuilt one by one"
else
  echo "FAIL twopass: $(log twopass) should report new junctions and lack: $BLOCK"; rc=1
fi
if [ "$THREADS" -gt 1 ]; then
  SHARED="SA entries shared among $THREADS threads in"
  for cfg in encode twopass; do
    if grep -qF -- "$SHARED" "$(log $cfg)" 2>/dev/null; then echo "PASS $cfg log: $SHARED"; else echo "FAIL $cfg: $(log $cfg) lacks: $SHARED"; rc=1; fi
  done
fi
exit $rc
