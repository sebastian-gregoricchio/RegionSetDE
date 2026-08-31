# .flattenRegionSets

Turns a collection of region sets into a single `GRanges`, one element
per row of the future counts matrix. The set name and the region
identifier are stored in the metadata columns, and the regions are
optionally cut into tiles.

## Usage

``` r
.flattenRegionSets(
  regionSet,
  tileWidth = NULL,
  partialTiles = TRUE,
  verbose = TRUE
)
```

## Arguments

- regionSet:

  `RegionSetDE` object, `GRangesList` or named list of `GRanges`.

- tileWidth:

  Numeric value with the width of the tiles. Default: `NULL`, one row
  per region.

- partialTiles:

  Logical value indicating whether the trailing shorter tile must be
  kept. Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A named `GRanges` with the `region.set`, `region.id` and `tile.id`
metadata columns.

## Author

Sebastian Gregoricchio
