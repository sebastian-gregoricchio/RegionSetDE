# .topVariableRows

Picks the most variable rows of a panel, which is what makes an
ordination read the structure between samples rather than the
differences in depth.

## Usage

``` r
.topVariableRows(counts, rowIndex, topRegions = 2000)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- rowIndex:

  Integer vector with the rows of the panel.

- topRegions:

  Numeric value with the number of rows kept, or `NULL` for all of them.

## Value

An integer vector with the rows kept.

## Author

Sebastian Gregoricchio
