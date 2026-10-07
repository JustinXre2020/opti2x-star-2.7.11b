#!/usr/bin/env bash
# End-to-end oracle: the candidate STAR must produce byte-identical outputs (SAM body, transcriptome BAM body,
# chimeric junctions/SAM, SJ.out.tab, ReadsPerGene, mapping stats) to upstream STAR 2.7.11b for every config.
# Usage: STAR_IDX=... STAR_R1=... STAR_R2=... STAR_GTF=... oracle_test.sh <candidate STAR> <upstream STAR> <workdir> [threads] [config ...]
# Default configs: default encode se_short chim. Upstream results are cached in <workdir>/upstream-<config>-t<threads>.
# STAR_CAND_ARGS, if set, is appended to the candidate's command line only.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
CAND=$1; UPS=$2; WORK=$3; THREADS=${4:-1}
shift $(( $# < 4 ? $# : 4 ))
CONFIGS=("$@")
[ ${#CONFIGS[@]} -gt 0 ] || CONFIGS=(default encode se_short chim)
mkdir -p "$WORK"

rc=0
for cfg in "${CONFIGS[@]}"; do
  ref=$WORK/upstream-$cfg-t$THREADS
  key="$(cat "$UPS" "$HERE/star_run.sh" | shasum -a 256 | cut -c1-16) $STAR_IDX $STAR_R1 ${STAR_R2:-} ${STAR_GTF:-}"
  if [ ! -f "$ref/oracle.txt" ] || [ "$(cat "$ref/cache_key" 2>/dev/null)" != "$key" ]; then
    if ! STAR_EXTRA_ARGS= "$HERE/star_run.sh" "$UPS" "$ref" "$cfg" "$THREADS"; then
      echo "FAIL $cfg: upstream run failed, see $ref/Log.out"; rc=1; continue
    fi
    echo "$key" > "$ref/cache_key"
  fi
  cand=$WORK/candidate-$cfg-t$THREADS
  if ! STAR_EXTRA_ARGS="${STAR_CAND_ARGS:-}" "$HERE/star_run.sh" "$CAND" "$cand" "$cfg" "$THREADS"; then
    echo "FAIL $cfg: candidate run failed, see $cand/Log.out"; rc=1; continue
  fi
  if diff "$ref/oracle.txt" "$cand/oracle.txt" > "$cand/oracle.diff"; then
    echo "PASS $cfg (upstream $(cut -d' ' -f1 "$ref/time.txt")s, candidate $(cut -d' ' -f1 "$cand/time.txt")s)"
  else
    echo "FAIL $cfg"; sed 's/^/    /' "$cand/oracle.diff"; rc=1
  fi
done
exit $rc
