# .storedRegionIds

Returns the identifier a counted object already carries for every
region, the `region.id` metadata column when the rows have one and their
names otherwise. Not to be confused with `.regionIdentifiers`, which
assigns the identifiers in the first place while a region set is loaded.

## Usage

``` r
.storedRegionIds(regions)
```

## Arguments

- regions:

  `GRanges` with the regions of the object.

## Value

A character vector with one identifier per region.

## Author

Sebastian Gregoricchio
