# .loadRegionList

Picks one of the region lists shipped with the package and reads it into
a `GRanges`.

## Usage

``` r
.loadRegionList(
  type,
  genome,
  assay = NULL,
  source = NULL,
  seqlevelsStyle = "UCSC",
  verbose = TRUE
)
```

## Arguments

- type:

  String with the type of list, `"blacklist"` or `"greenlist"`.

- genome:

  String with the genome assembly.

- assay:

  String with the assay the list was built for, or `NULL` for the lists
  tied to none.

- source:

  String with the source of the list, or `NULL` when the genome and the
  assay are enough to tell them apart.

- seqlevelsStyle:

  String with the chromosome naming style of the output, or `NULL`.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A `GRanges` with the regions of the list.

## Author

Sebastian Gregoricchio
