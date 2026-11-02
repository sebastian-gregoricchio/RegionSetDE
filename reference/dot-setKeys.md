# .setKeys

Builds the key identifying a region row inside a set, the set name and
the region identifier pasted together. It is what
[`geneSetsToRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/geneSetsToRegionSets.md)
returns and what
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
resolves a supplied set against, and it survives the collapsing of the
tiles, which row positions do not.

## Usage

``` r
.setKeys(x)
```

## Arguments

- x:

  Data.frame or `GRanges` carrying the `region.set` and `region.id`
  columns.

## Value

A character vector with one key per row.

## Author

Sebastian Gregoricchio
