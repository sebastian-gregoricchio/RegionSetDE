# .signalBrackets

Places the brackets joining the groups compared in a signal plot, and
writes the fold change and the adjusted p-value of the set on each of
them.

## Usage

``` r
.signalBrackets(
  plotTable,
  colTable,
  setResults,
  comparisons = NULL,
  contrastGroups = NULL,
  groupOrder,
  valueColumn = "mean.log2FC",
  perFacet = FALSE
)
```

## Arguments

- plotTable:

  Data.frame with the values being drawn.

- colTable:

  Data.frame with the samples and their groups, ordered.

- setResults:

  `RegionSetDE.setResults` object.

- comparisons:

  List of character vectors of length two, or `NULL`.

- contrastGroups:

  Character vector of length two with the levels the contrast compares,
  or `NULL`.

- groupOrder:

  Character vector with the groups, ordered.

- valueColumn:

  String with the fold change written on the bracket.

- perFacet:

  Logical value indicating whether the height must be computed inside
  each panel.

## Value

A data.frame with one row per bracket, or `NULL` when no comparison can
be drawn.

## Author

Sebastian Gregoricchio
