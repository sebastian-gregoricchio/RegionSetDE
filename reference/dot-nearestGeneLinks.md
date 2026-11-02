# .nearestGeneLinks

Builds the links of the regions by proximity to a transcription start
site, reading the sites from a `TxDb` or an `EnsDb`. Regions whose
nearest site sits further than `maxDistance` are written into the table
with a missing gene rather than left out of it, so that an unassigned
region stays distinguishable from a region nobody looked at.

## Usage

``` r
.nearestGeneLinks(
  regions,
  txdb = NULL,
  geneRanges = NULL,
  maxDistance = NULL,
  flankDistance = NULL,
  tssRegion = c(-3000, 3000),
  verbose = TRUE
)
```

## Arguments

- regions:

  `GRanges` with the regions of the object.

- txdb:

  `TxDb` or `EnsDb` object with the gene models. Default: `NULL`, the
  sites come from `geneRanges`.

- geneRanges:

  `GRanges` with the genes, or their start sites, carrying a gene
  identifier. Default: `NULL`, the genes come from `txdb`.

- maxDistance:

  Numeric value with the distance beyond which no gene is assigned, in
  base pairs. Default: `NULL`, every region takes its nearest site.

- flankDistance:

  Numeric value: every gene whose start site falls within it is linked,
  rather than the nearest one alone. Default: `NULL`, the nearest site.

- tssRegion:

  Numeric vector of length two with the window around a start site,
  signed and read on the strand of the gene, inside which a link counts
  as a promoter one. Default: `c(-3000, 3000)`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A data.frame in the long link format, with a `link.type` column holding
`"promoter"` or `"distal"`.

## Author

Sebastian Gregoricchio
