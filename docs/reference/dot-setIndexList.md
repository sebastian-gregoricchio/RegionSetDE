# .setIndexList

Resolves the sets a set level test runs on, from either the `region.set`
column of the rows or a named list of region keys supplied by the user.
The two differ in one way that matters: the sets of the column are a
partition, one row belonging to one set, while supplied sets overlap
freely, which is what gene sets do.

## Usage

``` r
.setIndexList(regionStats, regionSets = NULL)
```

## Arguments

- regionStats:

  Data.frame with the per-region statistics, carrying `region.set` and
  `region.id`.

- regionSets:

  Character vector with the names of the sets to test, or a named list
  of region keys. Default: `NULL`, every set of the column.

## Value

A named list with the row positions of every set.

## Author

Sebastian Gregoricchio
