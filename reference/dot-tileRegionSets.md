# .tileRegionSets

Cuts each region into adjacent tiles of fixed width, propagating the
metadata columns of the parent region to all its tiles.

## Usage

``` r
.tileRegionSets(regions, tileWidth, partialTiles = TRUE, verbose = TRUE)
```

## Arguments

- regions:

  `GRanges` with the flattened region sets.

- tileWidth:

  Numeric value with the width of the tiles, in base pairs.

- partialTiles:

  Logical value: `TRUE` keeps the trailing tile even when shorter than
  `tileWidth`, `FALSE` discards it. Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `GRanges` with one element per tile and an extra `tile.id` metadata
column.

## Author

Sebastian Gregoricchio
