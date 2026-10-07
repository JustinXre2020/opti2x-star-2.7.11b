#!/usr/bin/env bash
# --genomeMmap Yes must leave every output byte-identical to upstream, and must map SA/SAindex only when they are used as
# stored on disk: the default config maps them, the encode config (sjdb insertion from a GTF) keeps a private copy.
# Usage: STAR_IDX=... STAR_R1=... STAR_R2=... STAR_GTF=... test_genomeMmap.sh <candidate STAR> <upstream STAR> <workdir> [threads]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
CAND=$1; UPS=$2; WORK=$3; THREADS=${4:-1}

rc=0
bad=$WORK/candidate-genomeMmap-bad
rm -rf "$bad"; mkdir -p "$bad"
if "$CAND" --genomeDir "$STAR_IDX" --readFilesIn "$STAR_R1" --outFileNamePrefix "$bad/" --genomeMmap Maybe > "$bad/stdout.txt" 2>&1 \
   || ! grep -qF -- "--genomeMmap=Maybe" "$bad/stdout.txt"; then
  echo "FAIL --genomeMmap Maybe is not rejected, see $bad/stdout.txt"; rc=1
else
  echo "PASS --genomeMmap Maybe is rejected"
fi
STAR_CAND_ARGS="--genomeMmap Yes" "$HERE/oracle_test.sh" "$CAND" "$UPS" "$WORK" "$THREADS" default encode || rc=1

expectLog() {
  local log=$WORK/candidate-$1-t$THREADS/Log.out
  if grep -qF -- "$2" "$log" 2>/dev/null; then echo "PASS $1 log: $2"; else echo "FAIL $1: $log lacks: $2"; rc=1; fi
}
expectLog default "SA and SAindex are mapped read-only from the genome directory"
expectLog encode "--genomeMmap Yes is ignored"
exit $rc
