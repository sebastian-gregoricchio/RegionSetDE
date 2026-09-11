# .locateRegion

Finds the rows of a counts object belonging to one region.

## Usage

``` r
.locateRegion(counts, region)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- region:

  String with the region identifier, written as `"set|id"` or as the
  identifier alone, or a `GRanges` of length one.

## Value

A data.frame with the rows of the region and their positions in the
object.

## Author

Sebastian Gregoricchio
