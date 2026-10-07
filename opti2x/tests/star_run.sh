#!/usr/bin/env bash
# Usage: STAR_IDX=<genomeDir> STAR_R1=<r1.fq> STAR_R2=<r2.fq> STAR_GTF=<gtf> star_run.sh <STAR binary> <outdir> <config> [threads]
# STAR_EXTRA_ARGS, if set, is split on spaces and appended to the STAR command line.
# Runs one STAR alignment config and writes oracle hashes to <outdir>/oracle.txt
# and wall/RSS to <outdir>/time.txt. With more than one thread the order of reads in the outputs varies from run to run,
# so the SAM, BAM and chimeric junction records are hashed sorted; the primary transcriptome alignment of a read is
# drawn from a per-thread RNG, so then only the number of primary records of each read counts.
set -euo pipefail
BIN=$1; OUT=$2; CFG=$3; THREADS=${4:-1}
HERE=$(cd "$(dirname "$0")" && pwd)
: "${STAR_IDX:?set STAR_IDX to a STAR genome directory}"
: "${STAR_R1:?set STAR_R1/STAR_R2 to the read FASTQ files}"
IDX=$STAR_IDX
R1=$STAR_R1
R2=${STAR_R2:-}
GTF=${STAR_GTF:-}
rm -rf "$OUT"; mkdir -p "$OUT"

COMMON=(--runThreadN "$THREADS" --genomeDir "$IDX" --genomeLoad NoSharedMemory
        --outFileNamePrefix "$OUT/" --outSAMtype SAM --outSAMmode Full)
case "$CFG" in
  default)   ARGS=(--readFilesIn "$R1" "$R2" --outSAMattributes NH HI AS nM NM MD jM jI MC --outSAMstrandField intronMotif) ;;
  encode)    ARGS=(--readFilesIn "$R1" "$R2" --outFilterType BySJout --outSAMattributes NH HI AS NM MD
                   --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04
                   --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000
                   --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --peOverlapNbasesMin 10
                   --quantMode GeneCounts TranscriptomeSAM --sjdbGTFfile "${GTF:?encode needs STAR_GTF}") ;;
  se_short)  ARGS=(--readFilesIn "$R1" --clip3pNbases 38 --outFilterMismatchNoverLmax 0.1
                   --outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5 --outSAMattributes NH HI AS nM NM MD jM jI MC) ;;
  chim)      ARGS=(--readFilesIn "$R1" "$R2" --chimSegmentMin 12 --chimJunctionOverhangMin 8 --chimOutType Junctions SeparateSAMold
                   --alignEndsType EndToEnd --outSAMattributes NH HI AS nM) ;;
  twopass)   ARGS=(--readFilesIn "$R1" "$R2" --twopassMode Basic --outSAMattributes NH HI AS nM NM MD jM jI MC
                   --outSAMstrandField intronMotif) ;;
  *) echo "unknown config $CFG" >&2; exit 2 ;;
esac
EXTRA=(${STAR_EXTRA_ARGS:-})

/usr/bin/time -l "$BIN" "${COMMON[@]}" "${ARGS[@]}" ${EXTRA[@]+"${EXTRA[@]}"} > "$OUT/stdout.txt" 2> "$OUT/time_raw.txt"
WALL=$(awk '/real/ {print $1; exit}' "$OUT/time_raw.txt")
RSS=$(awk '/maximum resident set size/ {print $1; exit}' "$OUT/time_raw.txt")
echo "$WALL $RSS" > "$OUT/time.txt"

ordered() { if [ "$THREADS" -gt 1 ]; then LC_ALL=C sort; else cat; fi; }
BAMSORT=(); [ "$THREADS" -gt 1 ] && BAMSORT=(--sorted --anyPrimary)
{
  for f in Aligned.out.sam Aligned.toTranscriptome.out.bam Chimeric.out.junction Chimeric.out.sam SJ.out.tab ReadsPerGene.out.tab; do
    [ -f "$OUT/$f" ] || continue
    case "$f" in
      *.sam) h=$(grep -v '^@PG\|^@CO' "$OUT/$f" | ordered | shasum -a 256 | cut -c1-16) ;;
      *.bam) h=$(python3 "$HERE/bam_body_hash.py" ${BAMSORT[@]+"${BAMSORT[@]}"} "$OUT/$f") ;;
      *.junction) h=$(grep -v '^#' "$OUT/$f" | ordered | shasum -a 256 | cut -c1-16) ;;
      *) h=$(shasum -a 256 < "$OUT/$f" | cut -c1-16) ;;
    esac
    echo "$f $h"
  done
  # mapping statistics (skip timing/speed lines)
  echo "Log.final.out $(grep -v 'Started\|Finished\|speed' "$OUT/Log.final.out" | shasum -a 256 | cut -c1-16)"
} > "$OUT/oracle.txt"
