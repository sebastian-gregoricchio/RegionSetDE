# splitLoadRegions

Imports a collection of genomic regions stored in a single object or
file and splits them into individual region sets according to the values
of one of its columns. The regions can be provided as a `GRanges`, a
data.frame or the path to a BED-like/tabular file. All the arguments
controlling the import are passed to
[`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md).

## Usage

``` r
splitLoadRegions(
  regions,
  splitBy = "name",
  header = FALSE,
  selectedSets = NULL,
  minRegionsPerSet = 1,
  maxSets = 100,
  keepSplitColumn = FALSE,
  keepMetadata = TRUE,
  sortRegions = TRUE,
  reduceRegions = FALSE,
  removeDuplicatedRegions = TRUE,
  duplicatedSets = "stop",
  seqlevelsStyle = "UCSC",
  genomeAssembly = NULL,
  outputFormat = "RegionSetDE",
  verbose = TRUE
)
```

## Arguments

- regions:

  A `GRanges` object, a data.frame with chromosome/start/end columns, or
  a string indicating the path to a BED-like or tabular file containing
  all the regions.

- splitBy:

  String with the name of the column/metadata field containing the
  region set labels, or numeric value with its position in the input.
  The chromosome column (`"seqnames"`) is accepted as well, to split the
  regions by chromosome. Default: `"name"`, the fourth column of a BED
  file.

- header:

  Logical value to indicate whether the file contains a header line.
  Used only when `regions` is a file path. Notice that headerless files
  are read as BED and their coordinates are converted from 0-based to
  1-based, while data.frames and `GRanges` are assumed to be 1-based
  already. Default: `FALSE`.

- selectedSets:

  Character vector with the labels of the sets to import, the others
  being discarded. Default: `NULL`, all the sets are imported.

- minRegionsPerSet:

  Numeric value indicating the minimum number of regions that a set must
  contain to be imported. Sets falling below this threshold are dropped
  with a warning. Default: `1`.

- maxSets:

  Numeric value indicating the maximum number of sets tolerated.
  Splitting on a continuous column, such as the score, would generate
  thousands of single-region sets, therefore the loading is interrupted
  above this threshold. Default: `100`.

- keepSplitColumn:

  Logical value to indicate whether the splitting column must be kept in
  the metadata. Since it is constant within each set, it is redundant
  with the set name. Default: `FALSE`.

- keepMetadata:

  Logical value to indicate whether the metadata columns must be kept.
  Default: `TRUE`.

- sortRegions:

  Logical value to indicate whether the regions must be sorted by
  coordinate. Default: `TRUE`.

- reduceRegions:

  Logical value to indicate whether the overlapping regions within the
  same set must be collapsed
  ([`IRanges::reduce`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html)).
  Notice that metadata are lost upon reduction. Default: `FALSE`.

- removeDuplicatedRegions:

  Logical value to indicate whether the regions with identical
  coordinates within the same set must be collapsed to a single entry.
  Default: `TRUE`.

- duplicatedSets:

  String indicating how to handle two or more sets containing exactly
  the same regions: `"stop"`, `"remove"` or `"keep"`. Default: `"stop"`.

- seqlevelsStyle:

  String indicating the chromosome naming style to apply, one among
  `"UCSC"` (chr1), `"Ensembl"` (1) or `"NCBI"`. Default: `"UCSC"`.

- genomeAssembly:

  String indicating the genome assembly to store in the ranges seqinfo,
  e.g. `"hg38"`. Default: `NULL`, no assembly is assigned.

- outputFormat:

  String indicating the class of the returned object, one among
  `"RegionSetDE"`, `"GRangesList"` or `"list"`. Default:
  `"RegionSetDE"`.

- verbose:

  Logical value to indicate whether the loading messages must be
  printed. Default: `TRUE`.

## Value

A `RegionSetDE` object or, depending on `outputFormat`, a named
`GRangesList` or a named list of `GRanges`.

## See also

[`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
regionTable <- loadExampleData("regions", verbose = FALSE)

regionRanges <- GenomicRanges::makeGRangesFromDataFrame(regionTable,
                                                        keep.extra.columns = TRUE)

regions <- splitLoadRegions(regionRanges,
                            splitBy = "setName",
                            genomeAssembly = "rn4",
                            verbose = FALSE)
regionSetNames(regions)
#> [1] "promoterNonCpG" "intergenic"     "geneBody"       "promoterCpG"   

if (FALSE) { # \dontrun{
regions <- splitLoadRegions("peaks/all_peaks_annotated.bed",
                            splitBy = "name",
                            genomeAssembly = "hg38")
} # }
```
