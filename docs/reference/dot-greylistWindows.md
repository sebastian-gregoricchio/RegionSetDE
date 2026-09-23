# .greylistWindows

Cuts the chromosomes into windows of fixed width overlapping by half,
the last ones trimmed at the end of each chromosome.

## Usage

``` r
.greylistWindows(chromosomeLengths, binSize)
```

## Arguments

- chromosomeLengths:

  Named numeric vector with the length of every chromosome.

- binSize:

  Integer with the width of the windows.

## Value

A `GRanges` with the windows and the chromosome lengths in its
`seqinfo`.

## Author

Sebastian Gregoricchio
