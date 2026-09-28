# .collapseTileMatrix

Averages the rows of a value matrix over the tiles of each region.

## Usage

``` r
.collapseTileMatrix(expressionMatrix, tileMap, pooling = "mean")
```

## Arguments

- expressionMatrix:

  Numeric matrix with one row per tile.

- tileMap:

  Integer vector giving the collapsed row of every tile, as returned by
  `.collapseTileStats`.

- pooling:

  String with how the tiles are pooled, either `"mean"`, their average,
  or `"stouffer"`, their sum divided by the square root of their number,
  which keeps z-scores on the scale of a single one. Default: `"mean"`.

## Value

A numeric matrix with one row per region.

## Author

Sebastian Gregoricchio
