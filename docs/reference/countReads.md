# countReads

Counts the reads of a group of BAM files over the regions of a
`RegionSetDE` object. Paired-end data are counted as fragments, while
single-end reads are extended to the fragment length before the overlap
is evaluated. The regions can be recentred on the summit of the signal,
as DiffBind does with the peaks of a consensus, or cut into tiles of
fixed width, in which case each tile becomes a row of the resulting
object. The input libraries of the samples, when there are any, are
counted over the same rows and stored beside the counts.

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
  inputFiles = NULL,
  countInput = TRUE,
  summits = NULL,
  summitSource = "reads",
  nThreads = 1,
  verbose = TRUE
)
```

## Arguments

- regionSet:

  `RegionSetDE` object returned by
  [`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md)
  or
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md),
  a named `GRangesList`, or a single `GRanges`, which is loaded by
  [`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md)
  as one set named `regions`.

- bamFiles:

  Character vector with the paths of the BAM files. Each file must be
  indexed, and all of them must be aligned to the same assembly, whose
  chromosomes they may name in different styles (`chr1` in some files
  and `1` in others). Default: `NULL`, taken from `sampleSheet`, or from
  the sample sheet a consensus was built from by
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- sampleSheet:

  Data.frame returned by
  [`loadSampleSheet`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadSampleSheet.md),
  or the path to a sample sheet, providing the BAM files, the input
  files, the sample names and the annotation in one go. Default: `NULL`.

- sampleNames:

  Character vector with the sample names. Default: `NULL`, the BAM file
  names are used.

- sampleMetadata:

  Data.frame with the sample annotation, stored in the `colData`. When
  it contains a `sample` column the rows are matched by name, otherwise
  they must follow the order of `bamFiles`. Default: `NULL`.

- tileWidth:

  Numeric value with the width of the tiles, in base pairs. Cannot be
  combined with `summits`. Default: `NULL`, one row per region.

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

  Length to which single-end reads are extended. Either a number applied
  to every single-end sample; one number per BAM file, matched to the
  sample names when the vector is named; the name of a column of the
  sample sheet or of `sampleMetadata` holding one value per sample, such
  as the fragment length computed by phantompeakqualtools; or `"auto"`,
  to estimate it from the reads with
  [`estimateFragmentLength`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md).
  Ignored for paired-end samples, whose fragments are rebuilt from the
  pairs. Default: `150`.

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

- inputFiles:

  Character vector with the path of the input BAM file of every sample,
  `NA` for a sample without input. One input can serve several samples
  and is counted once. Default: `NULL`, the `input` column of the sample
  sheet or of `sampleMetadata`, when there is one.

- countInput:

  Logical value to indicate whether the input files must be counted.
  Default: `TRUE`.

- summits:

  Numeric value with the half width of the regions recentred on their
  summit, which then span `2 * summits + 1` bp, as with the `summits`
  argument of
  [`DiffBind::dba.count`](https://rdrr.io/pkg/DiffBind/man/dba.count.html).
  `0` locates the summits and stores them without moving the regions.
  Default: `NULL`, the regions are counted as they are.

- summitSource:

  String with where the summits are taken from: `"reads"`, the highest
  point of the fragment pileup of the samples, or `"peaks"`, the summits
  written by the peak caller in the narrowPeak files of a consensus
  built by
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).
  Default: `"reads"`.

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
`region.set` column of the `rowData`, and the length the single-end
reads were extended to in `fragment.length` (`NA` for paired-end
samples). With inputs, the `input` assay holds for every sample the
counts of its input over the same rows (`NA` for a sample without one),
and the `colData` gains `input.id` and `input.library.size`. With
`summits`, the `rowData` gains `summit`, the position of the summit of
every region.

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
reads extended to their fragment length, so that both end up with one
count per sequenced fragment. Forcing a paired-end file through the
single-end path counts each mate on its own and nearly doubles its
values, while the opposite mistake finds no pair and returns a column of
zeros, which is why the layout is read from the files by default. The
resolved layout of each sample is stored in the `paired.end` column of
the `colData`.

Single-end libraries rarely share the same fragment length, and a read
extended too far spills into the neighbouring regions while one extended
too little misses the centre of its own. The length can come from the
pipeline that produced the files, as a column of the sample sheet, or
from `fragmentLength = "auto"`, which runs the strand cross-correlation
of
[`estimateFragmentLength`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md)
over the regions being counted.
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md)
and
[`countGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countGreenlist.md)
reuse the lengths chosen here, sample by sample.

The inputs are counted with the same filters as the samples, over the
same rows, and a single-end input is extended to the mean fragment
length of the samples it serves. Nothing downstream subtracts them: the
counts of a sample stay the reads of that sample, as the count models
need. The input assay is there to check the enrichment of the regions,
and to see whether a change between conditions also shows in the inputs,
which points to copy number rather than to binding.
[`libInfo`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/libInfo.md)
reports the inputs beside the libraries.

A consensus of peaks has regions of every width, and a wide region
collects more background than a narrow one around the same summit. With
`summits`, every region is replaced by a window of fixed width centred
on its summit, which is what DiffBind does by default. With
`summitSource = "reads"` the summit is the middle of the highest stretch
of fragment pileup in each sample, averaged over the samples with
weights proportional to the height of their pileup, scaled by their
depth, so that the samples carrying signal decide where it sits. With
`summitSource = "peaks"` it is the average of the summits the peak
caller wrote for the peaks overlapping the region, weighted by their
significance, which needs no BAM file but works only for narrowPeak
files. A region with neither reads nor peaks keeps its midpoint. The
regions keep their identifiers, so each window can be traced back to the
region it came from, and windows of neighbouring regions may overlap, as
in DiffBind.

Regions and BAM files do not need to share the same chromosome naming
style, and neither do the BAM files among themselves. The names of the
first file are the reference. The regions, `discardRegions` and
`excludeChromosomes` are brought to them chromosome by chromosome, for
the counting only, and every other file is read under its own names.
UCSC regions can then be counted on Ensembl alignments, or on a mix of
the two, and the object still comes back with the names of the input
sets. The same holds for the inputs. What the files cannot differ in is
the assembly: two files giving different lengths to the same chromosome
are refused. Contigs that some files lack under any name, as scaffolds
and decoys often do between two builds of one assembly, hold no read in
those files and enter the library sizes of the others, unless they are
listed in `excludeChromosomes`.

A single `GRanges` is taken as one set of regions. It goes through
[`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md)
with its order and its chromosome names kept, the regions with identical
coordinates collapsed into one, and the set is called `regions`. To give
the set another name, to sort it or to split it into several sets, call
[`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md)
or
[`splitLoadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/splitLoadRegions.md)
first.

## See also

[`estimateFragmentLength`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md),
[`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md),
[`loadCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadCounts.md),
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md),
[`libInfo`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/libInfo.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
# The peaks of one sample of the AR example stand in for a region set
sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)

# The sheet brings the BAM files, the inputs, the sample names and the annotation
counts <- countReads(peakRegions, sampleSheet = sampleSheet, verbose = FALSE)
counts
#> class: RegionSetDE.counts 
#> dim: 101 9 
#> metadata(2): signal.type count.like
#> assays(2): counts input
#> rownames(101): peaks|chr19:46089737-46089961
#>   peaks|chr19:46300828-46301307 ... peaks|chr19:57823680-57823967
#>   peaks|chr19:57840110-57840669
#> rowData names(9): region.set region.id ... V9 V10
#> colnames(9): AR_DMSO_r1 AR_DMSO_r2 ... AR_R1881_24h_r2 AR_R1881_24h_r3
#> colData names(13): sample bam.file ... library.size input.library.size
head(SummarizedExperiment::assay(counts, "input"), 3)
#>                               AR_DMSO_r1 AR_DMSO_r2 AR_DMSO_r3 AR_R1881_4h_r1
#> peaks|chr19:46089737-46089961          0          0          0              0
#> peaks|chr19:46300828-46301307          1          1          1              1
#> peaks|chr19:46314001-46314455          0          0          0              0
#>                               AR_R1881_4h_r2 AR_R1881_4h_r3 AR_R1881_24h_r1
#> peaks|chr19:46089737-46089961              0              0               0
#> peaks|chr19:46300828-46301307              1              1               1
#> peaks|chr19:46314001-46314455              0              0               0
#>                               AR_R1881_24h_r2 AR_R1881_24h_r3
#> peaks|chr19:46089737-46089961               0               0
#> peaks|chr19:46300828-46301307               1               1
#> peaks|chr19:46314001-46314455               0               0

# The same files given one by one, with the annotation as a table and no input
counts <- countReads(peakRegions,
                     bamFiles = sampleSheet$bam,
                     sampleNames = sampleSheet$sample,
                     sampleMetadata = sampleSheet[, c("sample", "condition")],
                     verbose = FALSE)

# Windows of 401 bp centred on the summit of the reads, as DiffBind counts a consensus
summitCounts <- countReads(peakRegions, sampleSheet = sampleSheet, summits = 200, verbose = FALSE)
head(SummarizedExperiment::rowRanges(summitCounts), 3)
#> GRanges object with 3 ranges and 10 metadata columns:
#>                                 seqnames            ranges strand |  region.set
#>                                    <Rle>         <IRanges>  <Rle> | <character>
#>   peaks|chr19:46089737-46089961    chr19 46089699-46090099      * |       peaks
#>   peaks|chr19:46300828-46301307    chr19 46300862-46301262      * |       peaks
#>   peaks|chr19:46314001-46314455    chr19 46313939-46314339      * |       peaks
#>                                              region.id   tile.id
#>                                            <character> <integer>
#>   peaks|chr19:46089737-46089961 chr19:46089737-46089..      <NA>
#>   peaks|chr19:46300828-46301307 chr19:46300828-46301..      <NA>
#>   peaks|chr19:46314001-46314455 chr19:46314001-46314..      <NA>
#>                                                   name     score        V7
#>                                            <character> <integer> <numeric>
#>   peaks|chr19:46089737-46089961 AR_R1881_24h_r1_peak..       159   9.71586
#>   peaks|chr19:46300828-46301307 AR_R1881_24h_r1_peak..      1963  46.75880
#>   peaks|chr19:46314001-46314455 AR_R1881_24h_r1_peak..       149   8.12502
#>                                        V8        V9       V10    summit
#>                                 <numeric> <numeric> <integer> <integer>
#>   peaks|chr19:46089737-46089961   18.4847   15.9682       103  46089899
#>   peaks|chr19:46300828-46301307  199.8470  196.3410       236  46301062
#>   peaks|chr19:46314001-46314455   17.4005   14.9008       186  46314139
#>   -------
#>   seqinfo: 1 sequence from hg38 genome; no seqlengths

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
