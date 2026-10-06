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

  `GRanges` with the windows, named as in `chromosomeLengths`.

- chromosomeLengths:

  Named numeric vector with the length of every chromosome, as returned
  by `.bamChromosomeMap`.

- ...:

  Read filters passed to `.windowFragments`.

## Value

A list with `coverage`, an `RleList` with one element per chromosome of
`chromosomeLengths`, and `fragments`, the number of fragments read.

## Author

Sebastian Gregoricchio
