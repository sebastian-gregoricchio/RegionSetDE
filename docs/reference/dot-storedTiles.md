# .storedTiles

Internal accessor of the tile level table of a per-region result. Unlike
[`tileTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/tileTable.md)
it does not stop when the counts were not tiled.

## Usage

``` r
.storedTiles(results)
```

## Arguments

- results:

  `RegionSetDE.results` object.

## Value

A data.frame, empty when the counts were not tiled.

## Author

Sebastian Gregoricchio
