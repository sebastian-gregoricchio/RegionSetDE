# .carryLinkColumns

Copies the representative link columns from the annotated rows onto
another set of ranges describing the same regions, matched on the
coordinates.

## Usage

``` r
.carryLinkColumns(regions, annotatedRegions)
```

## Arguments

- regions:

  `GRanges` to write the columns onto.

- annotatedRegions:

  `GRanges` carrying the link columns.

## Value

The input ranges with the link columns added.

## Author

Sebastian Gregoricchio
