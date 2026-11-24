#!/usr/bin/env bash

## CRAM example of RegionSetDE, read by the tests of the CRAM reader.
##
## The alignments are those of ex1.bam, the paired-end example of Rsamtools (two contigs of the samtools
## examples), converted to CRAM against ex1.fa, the reference of the samtools sources
## (https://github.com/samtools/samtools/tree/develop/examples). Both contigs are written in one
## container (multi_seq_per_slice=1): this is the layout in which htslib needs the reference to rebuild
## the template length, so the tests reading the file with no reference at hand cover the hardest case.
## The reference paths (UR) and the program lines (@PG) are dropped from the header, which keeps only
## the names, the lengths and the MD5 sums of the contigs.
##
## overrunPairs.sam holds six pairs written by hand on a 2 kb contig whose sequence repeats ACGTTGCA,
## one per case of the fragment rule: a normal pair, a read running one base past the 5' end of its mate,
## two mates run past each other, a complete overrun, two mates on the same strand and a first mate on
## the reverse strand. The tests build its BAM file with Rsamtools::asBam() and read overrunPairs.cram,
## converted below, with no reference at hand.
##
## Made with samtools 1.19.2 (htslib 1.19). Run it from the root of the package.

set -euo pipefail

BAM=$(Rscript -e 'cat(system.file("extdata", "ex1.bam", package = "Rsamtools"))')
OUT=inst/extdata/cram
TMP=$(mktemp -d)

mkdir -p $OUT

#---------------------------------------#
# Reference of the example              #
#---------------------------------------#
curl -sL https://raw.githubusercontent.com/samtools/samtools/develop/examples/ex1.fa -o $TMP/ex1.fa
samtools faidx $TMP/ex1.fa


#---------------------------------------#
# Conversion and clean header           #
#---------------------------------------#
samtools view --no-PG -C -T $TMP/ex1.fa --output-fmt-option multi_seq_per_slice=1 -o $TMP/ex1.cram $BAM
samtools view -H --no-PG $TMP/ex1.cram | sed 's/\tUR:[^\t]*//' > $TMP/header.sam
samtools reheader -P $TMP/header.sam $TMP/ex1.cram > $OUT/ex1.cram
samtools index $OUT/ex1.cram



#---------------------------------------#
# Designed pairs                        #
#---------------------------------------#
awk 'BEGIN{print ">chrT"; s=""; for(i=0;i<250;i++) s=s "ACGTTGCA"; for(i=1;i<=length(s);i+=70) print substr(s,i,70)}' > $TMP/chrT.fa
samtools faidx $TMP/chrT.fa
samtools view -b -o $TMP/overrunPairs.bam $OUT/overrunPairs.sam
samtools view --no-PG -C -T $TMP/chrT.fa -o $TMP/overrunPairs.cram $TMP/overrunPairs.bam
samtools view -H --no-PG $TMP/overrunPairs.cram | sed 's/\tUR:[^\t]*//' > $TMP/overrunHeader.sam
samtools reheader -P $TMP/overrunHeader.sam $TMP/overrunPairs.cram > $OUT/overrunPairs.cram
samtools index $OUT/overrunPairs.cram

rm -r $TMP
