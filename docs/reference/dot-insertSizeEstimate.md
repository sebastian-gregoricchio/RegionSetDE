# .insertSizeEstimate

Reports the median insert size of the proper pairs of one paired-end BAM
file over a set of regions.

## Usage

``` r
.insertSizeEstimate(
  bamFile,
  regions,
  chromosomeLengths,
  minMapq,
  removeDuplicates,
  discardRegions
)
```

## Arguments

- bamFile:

  String with the path of the BAM file.

- regions:

  `GRanges` with the regions, named after the chromosomes of the BAM
  file.

- chromosomeLengths:

  Named numeric vector with the length of every chromosome of the BAM
  file.

- minMapq:

  Numeric value with the minimum mapping quality of a read.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded.

- discardRegions:

  `GRanges` with the regions whose reads must be ignored, or `NULL`.

## Value

A list with `read.length` (`NA`), `fragment.length` and `n.reads`, the
number of pairs.

## Author

Sebastian Gregoricchio
