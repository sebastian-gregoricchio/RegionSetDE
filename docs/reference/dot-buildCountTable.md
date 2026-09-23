# .buildCountTable

Builds the table returned by `countTable` from a `RegionSetDE.counts`
object.

## Usage

``` r
.buildCountTable(
  counts,
  level = "region",
  normalized = FALSE,
  format = "wide",
  set = NULL,
  tileSummary = NULL,
  extraColumns = TRUE,
  verbose = TRUE
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- level:

  String, either `"region"` or `"tile"`.

- normalized:

  Logical value indicating whether the normalised assay must be read.

- format:

  String, one among `"wide"`, `"long"` and `"matrix"`.

- set:

  Character vector with the region sets to keep, or `NULL`.

- tileSummary:

  String with the rule combining the tiles, or `NULL`.

- extraColumns:

  `TRUE`, `FALSE` or a character vector with the annotation columns
  wanted.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A data.frame or a numeric matrix, see `countTable`.

## Author

Sebastian Gregoricchio
