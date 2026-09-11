# .panelSets

Splits the rows of a counts object into the panels a sample level plot
will draw.

## Usage

``` r
.panelSets(counts, set = NULL, facetBySet = FALSE)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- set:

  Character vector with the region sets to keep, or `NULL`.

- facetBySet:

  Logical value indicating whether every set gets a panel of its own.

## Value

A named list of integer vectors, one per panel.

## Author

Sebastian Gregoricchio
