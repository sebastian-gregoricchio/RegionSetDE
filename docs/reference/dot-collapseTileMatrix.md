# .collapseTileMatrix

Averages the rows of a value matrix over the tiles of each region.

## Usage

``` r
.collapseTileMatrix(expressionMatrix, tileMap)
```

## Arguments

- expressionMatrix:

  Numeric matrix with one row per tile.

- tileMap:

  Integer vector giving the collapsed row of every tile, as returned by
  `.collapseTileStats`.

## Value

A numeric matrix with one row per region.

## Author

Sebastian Gregoricchio
