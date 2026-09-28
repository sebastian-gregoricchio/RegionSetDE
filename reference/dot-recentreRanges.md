# .recentreRanges

Builds windows of fixed width around the summits, clipped to the
chromosome ends.

## Usage

``` r
.recentreRanges(summitPosition, halfWidth, chromosomeLength = NULL)
```

## Arguments

- summitPosition:

  Integer vector with the summits.

- halfWidth:

  Numeric value with the number of base pairs on each side of the
  summit.

- chromosomeLength:

  Numeric vector with the length of the chromosome of every summit, `NA`
  when unknown.

## Value

An `IRanges` with one window per summit.

## Author

Sebastian Gregoricchio
