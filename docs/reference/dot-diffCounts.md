# .diffCounts

Counts the changing regions of every set, for the annotation written in
the corners of the panels.

## Usage

``` r
.diffCounts(regionTable, bySet = TRUE)
```

## Arguments

- regionTable:

  Data.frame with the `region.set` and `diff.status` columns.

- bySet:

  Logical value indicating whether the counts must be split by region
  set.

## Value

A data.frame with the `n.up` and `n.down` columns, carrying `region.set`
when `bySet` is `TRUE`.

## Author

Sebastian Gregoricchio
