# .greenlistMatrix

Returns the matrix of greenlist counts the factors are estimated from,
either the one stored by `countGreenlist` or the one handed to
`normalizeCounts`.

## Usage

``` r
.greenlistMatrix(counts, greenlistCounts = NULL)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- greenlistCounts:

  Numeric matrix with one row per region and one column per sample, a
  numeric vector with one total per sample, or `NULL`.

## Value

A numeric matrix with one column per sample.

## Author

Sebastian Gregoricchio
