# countBigwig

Summarises the signal of a group of bigWig files over the regions of a
`RegionSetDE` object. Useful when the BAM files are not available, or
when the coverage has been produced by an external pipeline. The regions
can be cut into tiles of fixed width, in which case each tile becomes a
row of the resulting object.

## Usage

``` r
countBigwig(
  regionSet,
  bigwigFiles,
  sampleNames = NULL,
  sampleMetadata = NULL,
  tileWidth = NULL,
  partialTiles = TRUE,
  summaryFunction = "sum",
  missingAsZero = TRUE,
  roundValues = TRUE,
  nThreads = 1,
  verbose = TRUE
)
```

## Arguments

- regionSet:

  `RegionSetDE` object returned by
  [`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md),
  or a named `GRangesList`.

- bigwigFiles:

  Character vector with the paths of the bigWig files.

- sampleNames:

  Character vector with the sample names. Default: `NULL`, the bigWig
  file names are used.

- sampleMetadata:

  Data.frame with the sample annotation, stored in the `colData`. When
  it contains a `sample` column the rows are matched by name, otherwise
  they must follow the order of `bigwigFiles`. Default: `NULL`.

- tileWidth:

  Numeric value with the width of the tiles, in base pairs. Default:
  `NULL`, one row per region.

- partialTiles:

  Logical value: `TRUE` keeps the trailing tile of each region even when
  narrower than `tileWidth`, `FALSE` discards it together with the
  regions narrower than a single tile. Default: `TRUE`.

- summaryFunction:

  String indicating how the per-base values are collapsed into a single
  value per region, one among `"sum"`, `"mean"`, `"max"` or `"min"`.
  Default: `"sum"`.

- missingAsZero:

  Logical value indicating whether the positions not covered by the
  bigWig must be treated as zeros rather than as missing values.
  Default: `TRUE`.

- roundValues:

  Logical value indicating whether the summarised values must be rounded
  to integers. Default: `TRUE`.

- nThreads:

  Number of threads used to process the files in parallel. Default: `1`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.counts` object with one row per region, or per tile, and
one column per sample.

## Details

A bigWig holds coverage, not reads, so the library sizes cannot be
recovered from it: the `library.size` column of the `colData` is left as
`NA` and the total signal falling in the regions is reported in
`total.signal` instead. Normalisation factors must therefore be supplied
externally, or estimated from a background bigWig, and the values are
rounded by default because the count-based models expect integers. Files
carrying an already normalised coverage will produce values that no
longer follow a count distribution, which is worth keeping in mind at
the testing step.

## See also

[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
[`loadCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadCounts.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
# rtracklayer ships a small bigWig, and the regions are taken from its own content
bigwigFile <- file.path(system.file("tests", package = "rtracklayer"), "test.bw")

# The UCSC library behind rtracklayer reads a Windows drive letter as a URL
# protocol, so an absolute path to the packaged file cannot be opened there
if (.Platform$OS.type != "windows") {

  exampleRegions <- GenomicRanges::reduce(rtracklayer::import(bigwigFile))
  exampleRegions$setName <- "covered"

  exampleSets <- splitLoadRegions(exampleRegions, splitBy = "setName",
                                  seqlevelsStyle = NULL, verbose = FALSE)

  signal <- countBigwig(exampleSets,
                        bigwigFiles = bigwigFile,
                        sampleNames = "example",
                        verbose = FALSE)

  print(signal)
  print(SummarizedExperiment::assay(signal, "counts"))
}
#> class: RegionSetDE.counts 
#> dim: 2 1 
#> metadata(1): signal.type
#> assays(1): counts
#> rownames(2): covered|chr2:1-1500 covered|chr19:1501-2700
#> rowData names(3): region.set region.id tile.id
#> colnames(1): example
#> colData names(4): sample bigwig.file library.size total.signal
#>                         example
#> covered|chr2:1-1500        -750
#> covered|chr19:1501-2700     750
```
