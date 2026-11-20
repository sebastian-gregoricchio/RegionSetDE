# .windowCoverage

Reads the fragments of one BAM file around a set of windows and returns
their coverage, chromosome by chromosome.

## Usage

``` r
.windowCoverage(
  bamFile,
  windows,
  chromosomeLengths,
  isPairedEnd,
  fragmentLength,
  maxFragmentLength,
  minMapq,
  removeDuplicates,
  discardRegions = NULL
)
```

## Arguments

- bamFile:

  String with the path of the BAM file.

- windows:

  `GRanges` with the windows, named as in `chromosomeLengths`.

- chromosomeLengths:

  Named numeric vector with the length of every chromosome, as returned
  by `.bamChromosomeMap`.

- isPairedEnd:

  Logical value indicating whether the file holds paired-end reads.

- fragmentLength:

  Numeric value with the length single-end reads are extended to.

- maxFragmentLength:

  Numeric value with the maximum length of a paired-end fragment.

- minMapq:

  Numeric value with the minimum mapping quality of a read.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded.

- discardRegions:

  `GRanges` with the regions whose reads must be ignored, named as in
  `chromosomeLengths`. Default: `NULL`.

## Value

A list with `coverage`, an `RleList` with one element per chromosome of
`chromosomeLengths`, and `fragments`, the number of fragments read.

## Author

Sebastian Gregoricchio
