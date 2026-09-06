# .statsRanges

Rebuilds a `GRanges` from the coordinate columns of a statistics table.

## Usage

``` r
.statsRanges(regionStats, index = NULL)
```

## Arguments

- regionStats:

  Data.frame returned by `.setStatistics`.

- index:

  Integer vector with the rows to take. Default: `NULL`, all of them.

## Value

A `GRanges`.

## Author

Sebastian Gregoricchio
