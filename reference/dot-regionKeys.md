# .regionKeys

Builds the coordinate key used to match a link table to the rows of an
object. Coordinates rather than identifiers, because a link table
computed elsewhere knows where a region is and not what the package
decided to call it. The same coordinates appearing in two sets return
the same key twice, so one link reaches both rows, which is what a
region shared between two sets should do.

## Usage

``` r
.regionKeys(x)
```

## Arguments

- x:

  `GRanges` with the regions, or a data.frame carrying `seqnames`,
  `start` and `end`.

## Value

A character vector with one key per region, written as `chr:start-end`.

## Author

Sebastian Gregoricchio
