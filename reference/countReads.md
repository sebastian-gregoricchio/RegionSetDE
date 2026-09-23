# countReads

Counts the reads of a group of BAM files over the regions of a
`RegionSetDE` object. Paired-end data are counted as fragments, while
single-end reads are extended to the expected fragment length before the
overlap is evaluated. The regions can be cut into tiles of fixed width,
in which case each tile becomes a row of the resulting object.

## Usage

``` r
countReads(
  regionSet,
  bamFiles = NULL,
  sampleSheet = NULL,
  sampleNames = NULL,
  sampleMetadata = NULL,
  tileWidth = NULL,
  keepMetadata = TRUE,
  regionId = NULL,
  partialTiles = TRUE,
  pairedEnd = "auto",
  fragmentLength = 150,
  maxFragmentLength = 1000,
  minMapq = 20,
  removeDuplicates = TRUE,
  excludeChromosomes = NULL,
  discardRegions = NULL,
  fullLibrarySize = TRUE,
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
  indexed, and all of them must share the same header. Default: `NULL`,
  taken from `sampleSheet`, or from the sample sheet a consensus was
  built from by
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- sampleSheet:

  Data.frame returned by
  [`loadSampleSheet`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadSampleSheet.md),
  or the path to a sample sheet, providing the BAM files, the sample
  names and the annotation in one go. Default: `NULL`.

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

- pairedEnd:

  Logical value, one logical value per BAM file, or the string `"auto"`
  to read the layout from the files themselves. Default: `"auto"`.

- fragmentLength:

  Numeric value with the length to which single-end reads are extended.
  Applied to the single-end samples only. Default: `150`.

- maxFragmentLength:

  Numeric value with the maximum length accepted for a paired-end
  fragment. Applied to the paired-end samples only. Default: `1000`.

- minMapq:

  Numeric value with the minimum mapping quality of a read. Default:
  `20`.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded. Default: `TRUE`.

- excludeChromosomes:

  Character vector with the chromosomes left out of the library sizes,
  written in either naming style, `chrM` and `MT` both reaching the
  mitochondrial genome of a file naming it either way, for instance the
  mitochondrial genome, chrY or the unplaced and alternative contigs. A
  name matching no chromosome once converted raises a warning, since it
  would leave that chromosome inside the library sizes without a word.
  The regions lying on them are still counted: to leave those out as
  well, filter the regions when loading them. The contigs can be
  collected from the BAM header, e.g.
  `grep("_|EBV", names(Rsamtools::scanBamHeader(bamFile)[[1]]$targets), value = TRUE)`.
  Default: `NULL`, every chromosome enters the library sizes.

- discardRegions:

  `GRanges` with regions whose reads must be ignored, for instance a
  blacklist. A fragment is dropped when one of its reads starts inside
  them. Default: `NULL`.

- fullLibrarySize:

  Logical value: `TRUE` reads every chromosome that is not excluded,
  even those without any region, so that the library sizes cover the
  whole library; `FALSE` reads only the chromosomes carrying regions,
  which is much faster for a few regions but leaves library sizes that
  must not be used for normalisation. Default: `TRUE`.

- nThreads:

  Number of threads. The files are cut into pieces of at most 50 Mb,
  shared among the threads, so even a single file benefits from several
  of them. Default: `1`.

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

A paired-end fragment is counted in every region it overlaps, including
the regions it spans with both reads outside them. The fragment is
rebuilt from the first mate of each proper pair, whose position and
template length (TLEN) give its start and width, so the two reads never
have to be matched in memory. The pairs therefore have to be flagged as
proper by the aligner, and those longer than `maxFragmentLength` are
dropped. The mapping quality of the second mate is read from the `MQ`
tag, which `samtools fixmate` and Picard write; on files without it only
the first mate is checked, and a message says so.

Counts and library sizes are kept apart. The counts only need the
chromosomes carrying regions, and every region is counted, wherever it
lies. The library size of a sample is the number of fragments that went
through the same filters as the counts, on every chromosome of the BAM
files except those in `excludeChromosomes`. Leaving out the
mitochondrial genome matters in ATAC-seq, where its share of the reads
changes from sample to sample. With `fullLibrarySize = FALSE` only the
chromosomes carrying regions are read, and the library sizes are
partial: the counts do not change, but the library sizes are not usable
for normalisation, and
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)
warns when a method relies on them.

Paired-end and single-end samples can be mixed in the same call,
paired-end libraries being counted as fragments and single-end ones as
reads extended to `fragmentLength`, so that both end up with one count
per sequenced fragment. Forcing a paired-end file through the single-end
path counts each mate on its own and nearly doubles its values, while
the opposite mistake finds no pair and returns a column of zeros, which
is why the layout is read from the files by default. The resolved layout
of each sample is stored in the `paired.end` column of the `colData`.

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
# The peaks of one sample of the AR example stand in for a region set
sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)

# The sheet brings the BAM files, the sample names and the annotation
counts <- countReads(peakRegions, sampleSheet = sampleSheet, verbose = FALSE)
counts
#> class: RegionSetDE.counts 
#> dim: 101 9 
#> metadata(2): signal.type count.like
#> assays(1): counts
#> rownames(101): peaks|chr19:46089737-46089961
#>   peaks|chr19:46300828-46301307 ... peaks|chr19:57823680-57823967
#>   peaks|chr19:57840110-57840669
#> rowData names(9): region.set region.id ... V9 V10
#> colnames(9): AR_DMSO_r1 AR_DMSO_r2 ... AR_R1881_24h_r2 AR_R1881_24h_r3
#> colData names(11): sample bam.file ... paired.end library.size

# The same files given one by one, with the annotation as a table
counts <- countReads(peakRegions,
                     bamFiles = sampleSheet$bam,
                     sampleNames = sampleSheet$sample,
                     sampleMetadata = sampleSheet[, c("sample", "condition")],
                     verbose = FALSE)

# Tiles of 100 bp, one row each
tiledCounts <- countReads(peakRegions, sampleSheet = sampleSheet, tileWidth = 100, verbose = FALSE)
head(SummarizedExperiment::rowData(tiledCounts), 3)
#> DataFrame with 3 rows and 9 columns
#>                                      region.set              region.id
#>                                     <character>            <character>
#> peaks|chr19:46089737-46089961|tile1       peaks chr19:46089737-46089..
#> peaks|chr19:46089737-46089961|tile2       peaks chr19:46089737-46089..
#> peaks|chr19:46089737-46089961|tile3       peaks chr19:46089737-46089..
#>                                       tile.id                   name     score
#>                                     <integer>            <character> <integer>
#> peaks|chr19:46089737-46089961|tile1         1 AR_R1881_24h_r1_peak..       159
#> peaks|chr19:46089737-46089961|tile2         2 AR_R1881_24h_r1_peak..       159
#> peaks|chr19:46089737-46089961|tile3         3 AR_R1881_24h_r1_peak..       159
#>                                            V7        V8        V9       V10
#>                                     <numeric> <numeric> <numeric> <integer>
#> peaks|chr19:46089737-46089961|tile1   9.71586   18.4847   15.9682       103
#> peaks|chr19:46089737-46089961|tile2   9.71586   18.4847   15.9682       103
#> peaks|chr19:46089737-46089961|tile3   9.71586   18.4847   15.9682       103
```
