#!/bin/bash

set -ex
set -o pipefail

[ -z "$BREAKSEQ_PRE" ] && echo "BREAKSEQ_PRE not specified" && exit 1
[ -z "$INBAM" ] && echo "INBAM not specified" && exit 1
[ -z "$OUTBAM" ] && echo "OUTBAM not specified" && exit 1
[ -z "$BWA_BIN" -o -z "$SAMTOOLS_BIN" -o -z "$BPLIB" ] && echo "BWA_BIN, SAMTOOLS_BIN and BPLIB need to be specified" && exit 1

breakseq_opt_vars="min_soft_clip min_soft_clip_mapq min_soft_clip_mate_mapq bad_map_max_soft_clip bad_map_min_mapq bad_map_min_nm bad_map_min_mate_mapq"
breakseq_opt_str=""
for var in $breakseq_opt_vars
do
  [ -n "${!var}" ] && breakseq_opt_str="$breakseq_opt_str --$var ${!var}"
done

rm -f $OUTBAM
outfq=${OUTBAM/.bam}.fq.gz

echo "*** Extracting reads from $INBAM to $outfq ***"
[ ! -z "$CHR" ] && breakseq_opt_str="$breakseq_opt_str --chromosome $CHR"
PYTHONPATH="$PYTHONPATH" $BREAKSEQ_PRE $breakseq_opt_str --bams $INBAM | gzip -v -1 > $outfq

uncompressed_size=`LC_ALL=C gzip -l $outfq | awk 'NR==2 {print $2}'`
if ((uncompressed_size != 0)); then
  echo "*** Aligning reads from $outfq to breakpoint library $BPLIB ***"
  outsai=${OUTBAM/.bam}.sai
  $BWA_BIN aln -t 4 -q 10 $BPLIB $outfq > $outsai
  $BWA_BIN samse $BPLIB $outsai $outfq | $SAMTOOLS_BIN view -1 -F 4 -bo $OUTBAM -S -
  rm -f $outsai 
else
  echo "*** $outfq is empty so no alignment will be done ***"
fi

rm -f $outfq
