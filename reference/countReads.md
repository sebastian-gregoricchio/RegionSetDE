# countReads

Counts the reads of a group of BAM files over the regions of a
`RegionSetDE` object. Single-end reads are extended to the expected
fragment length before the overlap is evaluated, while paired-end data
are counted as fragments. The regions can be cut into tiles of fixed
width, in which case each tile becomes a row of the resulting object.

## Usage

``` r
countReads(
  regionSet,
  bamFiles,
  sampleNames = NULL,
  sampleMetadata = NULL,
  tileWidth = NULL,
  partialTiles = TRUE,
  pairedEnd = "auto",
  fragmentLength = 150,
  maxFragmentLength = 1000,
  minMapq = 20,
  removeDuplicates = TRUE,
  restrictChromosomes = NULL,
  discardRegions = NULL,
  nThreads = 1,
  verbose = TRUE
)
```

## Arguments

- regionSet:

  `RegionSetDE` object returned by
  [`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md),
  or a named `GRangesList`.

- bamFiles:

  Character vector with the paths of the BAM files. Each file must be
  indexed.

- sampleNames:

  Character vector with the sample names. Default: `NULL`, the BAM file
  names are used.

- sampleMetadata:

  Data.frame with the sample annotation, stored in the `colData`. When
  it contains a `sample` column the rows are matched by name, otherwise
  they must follow the order of `bamFiles`. Default: `NULL`.

- tileWidth:

  Numeric value with the width of the tiles, in base pairs. Default:
  `NULL`, one row per region.

- partialTiles:

  Logical value: `TRUE` keeps the trailing tile of each region even when
  narrower than `tileWidth`, `FALSE` discards it together with the
  regions narrower than a single tile. Default: `TRUE`.

- pairedEnd:

  Logical value, one logical value per BAM file, or the string `"auto"`
  to read the layout from the files themselves. Default: `"auto"`.

- fragmentLength:

  Numeric value with the length to which single-end reads are extended.
  Applied to the single-end samples only. Default: `150`.

- maxFragmentLength:

  Numeric value with the maximum insert size accepted for a pair.
  Applied to the paired-end samples only. Default: `1000`.

- minMapq:

  Numeric value with the minimum mapping quality of a read. Default:
  `20`.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded. Default: `TRUE`.

- restrictChromosomes:

  Character vector with the chromosomes to read from the BAM files.
  Default: `NULL`, all of them.

- discardRegions:

  `GRanges` with regions whose reads must be ignored, for instance a
  blacklist. Default: `NULL`.

- nThreads:

  Number of threads used to process the files in parallel. Default: `1`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.counts` object with one row per region, or per tile, and
one column per sample. The library sizes are stored in the
`library.size` column of the `colData`, the set membership in the
`region.set` column of the `rowData`.

## Details

Regions shared by several sets are counted only once and the values are
then copied to every set they belong to, which keeps the running time
proportional to the number of distinct regions.

Paired-end and single-end samples can be mixed in the same call. Each
layout is counted in a separate pass, paired-end libraries as fragments
and single-end ones as reads extended to `fragmentLength`, so that both
end up with one count per sequenced fragment. Forcing a paired-end file
through the single-end path counts each mate on its own and nearly
doubles its values, while the opposite mistake keeps only the proper
pairs and returns a column of zeros, which is why the layout is read
from the files by default. The resolved layout of each sample is stored
in the `paired.end` column of the `colData`.

Regions and BAM files do not need to share the same chromosome naming
style. When no chromosome is shared, the regions are converted to the
style of the files for the counting only, so that UCSC regions can be
counted on Ensembl alignments and the object still comes back with the
names of the input sets.

## See also

[`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md),
[`loadCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadCounts.md),
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
# Rsamtools ships a small alignment file, enough to run the counting itself
bamFile <- system.file("extdata", "ex1.bam", package = "Rsamtools")

exampleRegions <- GenomicRanges::GRanges(
  seqnames = rep(c("seq1", "seq2"), each = 3),
  ranges = IRanges::IRanges(start = rep(c(1, 500, 1000), 2), width = 300))

exampleRegions$setName <- rep(c("firstSet", "secondSet"), each = 3)

exampleSets <- splitLoadRegions(exampleRegions, splitBy = "setName",
                                seqlevelsStyle = NULL, verbose = FALSE)

counts <- countReads(exampleSets,
                     bamFiles = bamFile,
                     sampleNames = "example",
                     verbose = FALSE)
counts
#> class: RegionSetDE.counts 
#> dim: 6 1 
#> metadata(1): signal.type
#> assays(1): counts
#> rownames(6): firstSet|seq1:1-300 firstSet|seq1:500-799 ...
#>   secondSet|seq2:500-799 secondSet|seq2:1000-1299
#> rowData names(3): region.set region.id tile.id
#> colnames(1): example
#> colData names(4): sample bam.file paired.end library.size

SummarizedExperiment::assay(counts, "counts")
#>                          example
#> firstSet|seq1:1-300          113
#> firstSet|seq1:500-799        260
#> firstSet|seq1:1000-1299      277
#> secondSet|seq2:1-300         177
#> secondSet|seq2:500-799       323
#> secondSet|seq2:1000-1299     327

if (FALSE) { # \dontrun{
counts <- countReads(regions,
                     bamFiles = list.files("bam", pattern = "\\.bam$", full.names = TRUE),
                     sampleMetadata = data.frame(sample = c("ctrl1", "ctrl2", "treat1", "treat2"),
                                                 condition = c("ctrl", "ctrl", "treat", "treat")),
                     pairedEnd = TRUE,
                     nThreads = 4)

countsTiled <- countReads(regions, bamFiles = bamPaths, tileWidth = 500)
} # }

```
