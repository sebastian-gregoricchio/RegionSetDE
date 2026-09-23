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
  bigwigFiles = NULL,
  sampleSheet = NULL,
  sampleNames = NULL,
  sampleMetadata = NULL,
  tileWidth = NULL,
  keepMetadata = TRUE,
  regionId = NULL,
  partialTiles = TRUE,
  summaryFunction = "sum",
  missingAsZero = TRUE,
  countLike = FALSE,
  roundValues = FALSE,
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

  Character vector with the paths of the bigWig files. Default: `NULL`,
  taken from `sampleSheet`, or from the sample sheet a consensus was
  built from by
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- sampleSheet:

  Data.frame returned by
  [`loadSampleSheet`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadSampleSheet.md),
  or the path to a sample sheet, providing the bigWig files, the sample
  names and the annotation in one go. Default: `NULL`.

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

- keepMetadata:

  Logical value to indicate whether the metadata columns carried by the
  regions must be kept in the `rowData`, harmonised across the sets.
  Default: `TRUE`.

- regionId:

  String with the name of a metadata column holding the region
  identifiers, for instance a gene name. It must hold a different value
  for every region of every set. Default: `NULL`, the names of the
  ranges, and their coordinates when they are unnamed.

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

- countLike:

  Logical value with which you assert that the bigWig holds count-like
  values, meaning raw, unnormalised coverage that a count model may
  legitimately be applied to. Default: `FALSE`.

- roundValues:

  Logical value indicating whether the summarised values must be rounded
  to integers. Rounding is a formatting step and does not turn coverage
  into counts. Default: `FALSE`.

- nThreads:

  Number of threads used to process the files in parallel. Default: `1`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.counts` object with one row per region, or per tile, and
one column per sample. The `signal.type` entry of the metadata is set to
`"bigwig"` and `count.like` to whatever was declared through
`countLike`, which is what
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
reads to decide which engines it will run.

## Details

A bigWig holds coverage, not reads, so the library sizes cannot be
recovered from it: the `library.size` column of the `colData` is left as
`NA` and the total signal falling in the regions is reported in
`total.signal` instead. Normalisation factors must therefore be supplied
externally, or estimated from a background bigWig.

Coverage is not a count and rounding it does not make it one. A value of
12.72 becomes 13 and looks like a count, but the negative binomial
likelihood that `edgeR` and `DESeq2` are built on describes the number
of fragments falling in an interval, and a rounded coverage value is not
that number. The gap is widest for files carrying an already normalised
signal, CPM, RPKM, RPGC, fold enrichment over input, where the values
have been divided by a factor that the count model then has no way of
knowing about. It does not close entirely even for raw coverage:
coverage summed over an interval weights every fragment by how many of
its bases fall inside, so it is over-dispersed relative to the fragment
count it stands in for.

The object therefore records where its values came from, and
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
refuses the count engines on it unless `countLike` says otherwise. The
engine to reach for on bigWig input is `"limma"`, which models the log2
signal directly and asks nothing of the values that they cannot supply.
Setting `countLike = TRUE` is an assertion about the files, not a
setting: it says these are raw, unnormalised coverage tracks and the
count model is close enough for the purpose, and it should be stated in
the methods when it is used.

## See also

[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
[`loadCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadCounts.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
if (FALSE) { # \dontrun{
counts <- countBigwig(regions,
                      bigwigFiles = list.files("bigwig", pattern = "\\.bw$", full.names = TRUE),
                      summaryFunction = "sum",
                      nThreads = 4)
} # }
```
