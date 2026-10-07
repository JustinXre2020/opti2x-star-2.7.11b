#!/usr/bin/env bash
# Generating a genome with junctions from a GTF inserts them into the suffix array (sjdbInsertJunctions), and every index
# file must be byte-identical to upstream's. genomeParameters.txt is compared without its command-line lines.
# Usage: STAR_FASTA=<genome.fa> STAR_GTF=<gtf> test_genomeGenerate.sh <candidate STAR> <upstream STAR> <workdir> [threads]
# Upstream's index is cached in <workdir>/genome-upstream.
set -euo pipefail
: "${STAR_FASTA:?set STAR_FASTA to a genome FASTA}"
: "${STAR_GTF:?set STAR_GTF to a GTF}"
CAND=$1; UPS=$2; WORK=$3; THREADS=${4:-4}
mkdir -p "$WORK"

generate() {
  rm -rf "$2"; mkdir -p "$2"
  "$1" --runMode genomeGenerate --runThreadN "$THREADS" --genomeDir "$2" --genomeFastaFiles "$STAR_FASTA" \
       --genomeSAindexNbases 12 --sjdbGTFfile "$STAR_GTF" --sjdbOverhang 75 --outFileNamePrefix "$2/" > "$2/stdout.txt"
}

ref=$WORK/genome-upstream
key="$(shasum -a 256 < "$UPS" | cut -c1-16) $STAR_FASTA $STAR_GTF $THREADS"
if [ ! -f "$ref/SA" ] || [ "$(cat "$ref/cache_key" 2>/dev/null)" != "$key" ]; then
  generate "$UPS" "$ref"
  echo "$key" > "$ref/cache_key"
fi
cand=$WORK/genome-candidate
generate "$CAND" "$cand"

rc=0
for f in "$ref"/*; do
  name=$(basename "$f")
  case "$name" in
    Log.out|stdout.txt|cache_key) continue ;;
    genomeParameters.txt) same() { cmp -s <(grep -v '^###' "$1") <(grep -v '^###' "$2"); } ;;
    *) same() { cmp -s "$1" "$2"; } ;;
  esac
  if same "$f" "$cand/$name"; then echo "PASS $name"; else echo "FAIL $name differs"; rc=1; fi
done
grep -q "SA entries unchanged by the insertion" "$cand/Log.out" && { echo "FAIL a new genome's junctions cannot leave its SA unchanged"; rc=1; }
exit $rc
