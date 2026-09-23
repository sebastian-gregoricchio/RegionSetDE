#!/usr/bin/env bash

## Example BAM and peak files of RegionSetDE, cut out of the AR ChIP-seq of GSE284522
## (LNCaP, DMSO / R1881 4 h / R1881 24 h, three replicates each, plus one input; hg38, paired-end).
##
## The libraries themselves are far too heavy to ship, so what ends up in inst/extdata is a window of
## chromosome 19 subsampled to a fraction of its depth, with the sequence, the base qualities, the
## aligner tags, the read names and every contig the window does not touch thrown away. Counting only
## needs the coordinates and the CIGAR, so the counts over the window are the counts of the real
## libraries scaled by FRAC, and nothing else of the alignment survives. The peaks are the calls of
## the full libraries, filtered to the same window: they are not recalled on the slice, so they keep
## the depth the caller actually had.
##
## Run it from the folder holding the sorted BAM files and the narrowPeak files.

set -euo pipefail


#---------------------------------------#
# Settings                              #
#---------------------------------------#
CHROM=19                # naming as it is written in the BAM header, no chr prefix here
START=46000000
END=58000000

FRAC=0.12               # share of the read pairs kept, the knob for the file size
SEED=42                 # fixed so that a rerun cuts the same reads
THREADS=10
OUT=slice_data

WIN=${CHROM}:${START}-${END}

FILES=$(ls *Peak* | sed 's/.filtered.BAMPE_peaks.narrowPeak//')
LIBS="$FILES LNCaP_NT_input"

mkdir -p $OUT


#---------------------------------------#
# Window out of each library            #
#---------------------------------------#
# Cut once and kept, so that FRAC can be retuned below without touching the big files again
for i in $LIBS
do
  if [ ! -f $OUT/${i}_slice.bam ]
  then
    samtools view -@$THREADS -b -o $OUT/${i}_slice.bam ${i}_mapq20_mdup_sorted.bam $WIN
    samtools index -@$THREADS $OUT/${i}_slice.bam
  fi
done


#---------------------------------------#
# The file that gets shipped            #
#---------------------------------------#
printf 'library\treads_slice\treads_example\tsize_bam\tsize_bai\n' > $OUT/slice_report.tsv

for i in $LIBS
do
  # -s keeps or drops a whole template, so no mate is ever orphaned, and -f 2 leaves the proper pairs
  # alone, which is the same set the BAMPE calls were made on
  # On samtools older than 1.14 the seed goes in the fraction: -s ${SEED}.${FRAC#0.}
  samtools view -@$THREADS -h --no-PG -f 2 -s $FRAC --subsample-seed $SEED \
    --remove-tag MD,NM,MC,MQ,RG,XS,SA $OUT/${i}_slice.bam \
  | awk -v chrom="$CHROM" '
      BEGIN {FS = OFS = "\t"}

      # Of the header only the sort order and the one contig the reads sit on are worth keeping. The
      # rest is the alt and decoy list of hg38, which the index then carries an empty entry for.
      /^@HD/ {print; next}
      /^@SQ/ {if ($2 == "SN:" chrom) print; next}
      /^@/   {next}

      # Names go to r1, r2, and so on, since an Illumina name is the longest field left once the
      # sequence and the qualities are gone. Both mates keep the same one.
      {if (!($1 in name)) name[$1] = "r" (++k); $1 = name[$1]; $10 = "*"; $11 = "*"; print}' \
  | samtools view -@$THREADS -b --no-PG -o $OUT/${i}_example.bam

  samtools index -@$THREADS $OUT/${i}_example.bam

  samtools quickcheck $OUT/${i}_example.bam || { echo "broken file for $i"; exit 1; }

  # Where the weight sits and how many reads are left, to decide where FRAC should go
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$i" \
    "$(samtools view -c -@$THREADS $OUT/${i}_slice.bam)" \
    "$(samtools view -c -@$THREADS $OUT/${i}_example.bam)" \
    "$(du -h $OUT/${i}_example.bam | cut -f1)" \
    "$(du -h $OUT/${i}_example.bam.bai | cut -f1)" >> $OUT/slice_report.tsv
done

column -t $OUT/slice_report.tsv
du -ch $OUT/*_example.bam* | tail -n 1


#---------------------------------------#
# Peaks over the same window            #
#---------------------------------------#
for i in $FILES
do
  awk -v OFS="\t" -v chrom="$CHROM" -v start="$START" -v end="$END" \
    '$1 == chrom && $2 >= start && $3 <= end' \
    ${i}.filtered.BAMPE_peaks.narrowPeak > $OUT/${i}_peaks.narrowPeak
done

wc -l $OUT/*_peaks.narrowPeak
