# .windowCoverage

Reads the fragments of one BAM file around a set of windows and returns
their coverage, chromosome by chromosome.

## Usage

``` r
.windowCoverage(bamFile, windows, chromosomeLengths, ...)
```

## Arguments

- bamFile:

  String with the path of the BAM file.

- windows:

  `GRanges` with the windows, named after the chromosomes of the BAM
  file.

- chromosomeLengths:

  Named numeric vector with the length of every chromosome of the BAM
  file.

- ...:

  Read filters passed to `.windowFragments`.

## Value

A list with `coverage`, an `RleList` with one element per chromosome of
the BAM file, and `fragments`, the number of fragments read.

## Author

Sebastian Gregoricchio
