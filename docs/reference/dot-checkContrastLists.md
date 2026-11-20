# .checkContrastLists

Checks the list slots shared by the result classes: the cut-offs in
`thresholds` and the two groups of the contrast in `contrast.groups`.
Used by the validity methods.

## Usage

``` r
.checkContrastLists(object)
```

## Arguments

- object:

  `RegionSetDE.results`, `RegionSetDE.setResults` or
  `RegionSetDE.geneResults` object.

## Value

A character vector with the problems found, empty when there is none.

## Author

Sebastian Gregoricchio
