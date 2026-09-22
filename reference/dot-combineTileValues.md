# .combineTileValues

Combines the rows of a value matrix over the tiles of each region.

## Usage

``` r
.combineTileValues(valueMatrix, regionIndex, tileWidths, summaryRule = "sum")
```

## Arguments

- valueMatrix:

  Numeric matrix with one row per tile.

- regionIndex:

  Integer vector giving the region of every tile, numbered in order of
  first appearance.

- tileWidths:

  Numeric vector with the width of every tile.

- summaryRule:

  String, one among `"sum"`, `"mean"`, `"max"` and `"min"`.

## Value

A numeric matrix with one row per region, in the order of `regionIndex`.

## Author

Sebastian Gregoricchio
