# loadCounts

Imports a count matrix computed outside R, for instance by
featureCounts, bedtools multicov or deeptools multiBigwigSummary, and
attaches it to the regions of a `RegionSetDE` object. The rows of the
matrix are matched to the regions either by coordinates or by
identifier.

## Usage

``` r
loadCounts(
  regionSet,
  counts,
  sampleNames = NULL,
  sampleMetadata = NULL,
  countColumns = NULL,
  coordinateColumns = NULL,
  idColumn = NULL,
  matchBy = "coordinates",
  startsAt = 1,
  tileWidth = NULL,
  partialTiles = TRUE,
  missingRegions = "stop",
  librarySizes = NULL,
  header = TRUE,
  verbose = TRUE
)
```

## Arguments

- regionSet:

  `RegionSetDE` object returned by
  [`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md),
  or a named `GRangesList`.

- counts:

  Matrix, data.frame or path to a tab separated file with the counts.
  Lines starting with `#`, such as the featureCounts header, are
  skipped.

- sampleNames:

  Character vector with the sample names. Default: `NULL`, the names of
  the count columns are used.

- sampleMetadata:

  Data.frame with the sample annotation, stored in the `colData`. When
  it contains a `sample` column the rows are matched by name, otherwise
  they must follow the order of the count columns. Default: `NULL`.

- countColumns:

  Character vector or numeric positions of the columns holding the
  counts. Default: `NULL`, every numeric column that is not a coordinate
  or a standard annotation column.

- coordinateColumns:

  Character vector or numeric positions of the three columns holding
  chromosome, start and end. Used only when `matchBy = "coordinates"`.
  Default: `NULL`, detected from the column names.

- idColumn:

  String or numeric position of the column holding the region
  identifiers. Used only when `matchBy = "id"`. Default: `NULL`, the row
  names are used.

- matchBy:

  String indicating how the rows of `counts` are assigned to the
  regions, either `"coordinates"` or `"id"`. Default: `"coordinates"`.

- startsAt:

  Numeric value with the coordinate system of the count table, `1` for
  1-based tables such as featureCounts, `0` for BED-like tables.
  Default: `1`.

- tileWidth:

  Numeric value with the width of the tiles, to be set when the external
  matrix has been computed on tiles rather than on whole regions.
  Default: `NULL`.

- partialTiles:

  Logical value indicating whether the trailing shorter tile of each
  region has been kept. Default: `TRUE`.

- missingRegions:

  String indicating what to do with the regions absent from the count
  table, one among `"stop"`, `"zero"` or `"drop"`. Default: `"stop"`.

- librarySizes:

  Numeric vector with the library size of each sample, in the same order
  as the count columns. Default: `NULL`, the column sums of the imported
  table are used.

- header:

  Logical value indicating whether the file carries a header line.
  Ignored when `counts` is not a file path. Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.counts` object with one row per region, or per tile, and
one column per sample.

## Details

The column sums of the imported table are a poor substitute for the real
library sizes, since they only cover the regions present in the file.
When the sequencing depth is known it should be passed through
`librarySizes`, otherwise the normalisation should rely on factors
estimated elsewhere. Rows of the count table that match no region are
ignored, which makes it safe to import a genome wide matrix and keep
only the sets of interest.

## See also

[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
[`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
regionRanges <- SummarizedExperiment::rowRanges(counts)

# A count table as it would come out of featureCounts or a coverage tool
countTable <- data.frame(
  seqnames = as.character(GenomicRanges::seqnames(regionRanges)),
  start = GenomicRanges::start(regionRanges),
  end = GenomicRanges::end(regionRanges),
  SummarizedExperiment::assay(counts, "counts"),
  check.names = FALSE)

regions <- splitLoadRegions(regionRanges, splitBy = "region.set",
                            genomeAssembly = "rn4", verbose = FALSE)

reloaded <- loadCounts(regions, counts = countTable, verbose = FALSE)
reloaded
#> class: RegionSetDE.counts 
#> dim: 3224 4 
#> metadata(1): signal.type
#> assays(1): counts
#> rownames(3224): promoterNonCpG|promoterNonCpG|region_00002
#>   promoterNonCpG|promoterNonCpG|region_00003 ...
#>   promoterCpG|promoterCpG|region_03797
#>   promoterCpG|promoterCpG|region_03798
#> rowData names(3): region.set region.id tile.id
#> colnames(4): lv-H3K4me3-BN-female-bio1-tech1
#>   lv-H3K4me3-BN-male-bio2-tech1 lv-H3K4me3-SHR-male-bio2-tech1
#>   lv-H3K4me3-SHR-male-bio3-tech1
#> colData names(3): sample source.column library.size
```
