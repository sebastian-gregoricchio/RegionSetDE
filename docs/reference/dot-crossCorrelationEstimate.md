# .crossCorrelationEstimate

Estimates the fragment length of one single-end BAM file from the
distances between the 5' ends of its forward and reverse reads.

## Usage

``` r
.crossCorrelationEstimate(
  bamFile,
  regions,
  chromosomeLengths,
  maxDistance,
  minMapq,
  removeDuplicates,
  discardRegions,
  maxReads
)
```

## Arguments

- bamFile:

  String with the path of the BAM file.

- regions:

  `GRanges` with the regions, named after the chromosomes of the BAM
  file and without overlaps.

- chromosomeLengths:

  Named numeric vector with the length of every chromosome of the BAM
  file.

- maxDistance:

  Integer value with the longest distance considered.

- minMapq:

  Numeric value with the minimum mapping quality of a read.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded.

- discardRegions:

  `GRanges` with the regions whose reads must be ignored, or `NULL`.

- maxReads:

  Numeric value with the maximum number of forward reads used.

## Value

A list with `read.length`, `fragment.length`, `n.reads`, the number of
forward reads used, and `profile`, a data.frame with the number of pairs
at every distance.

## Author

Sebastian Gregoricchio
