# .windowFragments

Reads the fragments of one BAM file that overlap a set of windows, with
the same filters as the counting. The windows are widened by the longest
fragment accepted, so that a fragment reaching a window from outside is
not lost, and merged, so that every record is read once.

## Usage

``` r
.windowFragments(
  bamFile,
  windows,
  isPairedEnd,
  fragmentLength,
  maxFragmentLength,
  minMapq,
  removeDuplicates,
  chromosomeLengths,
  discardRegions = NULL,
  padding = NULL,
  readsOnly = FALSE
)
```

## Arguments

- bamFile:

  String with the path of the BAM file.

- windows:

  `GRanges` with the windows, named after the chromosomes of the BAM
  file.

- isPairedEnd:

  Logical value, `TRUE` for a paired-end file.

- fragmentLength:

  Numeric value with the length to which single-end reads are extended.

- maxFragmentLength:

  Numeric value with the maximum length of a paired-end fragment.

- minMapq:

  Numeric value with the minimum mapping quality of a read.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded.

- chromosomeLengths:

  Named numeric vector with the length of every chromosome of the BAM
  file.

- discardRegions:

  `GRanges` with the regions whose reads must be ignored, named after
  the chromosomes of the BAM file. Default: `NULL`.

- padding:

  Numeric value with the number of base pairs added on both sides of the
  windows. Default: `NULL`, the longest fragment accepted.

- readsOnly:

  Logical value: `TRUE` returns the reads as they are, with their strand
  and 5' end, instead of the fragments. Only for single-end files.
  Default: `FALSE`.

## Value

A `GRanges` with one element per fragment, or per read when
`readsOnly = TRUE`, carrying the 5' end of the read in the `five.prime`
column.

## Author

Sebastian Gregoricchio
