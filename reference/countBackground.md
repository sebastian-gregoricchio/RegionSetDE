# countBackground

Counts the reads falling in large genome wide bins, which provide the
background distribution used to estimate composition-aware normalisation
factors. The bins overlapping the counted regions are removed by
default, so that the normalisation is not driven by the signal under
study. The result is stored in the metadata of the counts object, where
the normalisation step retrieves it.

## Usage

``` r
countBackground(
  counts,
  bamFiles = NULL,
  binSize = 10000,
  excludeRegions = TRUE,
  minCount = 1,
  excludeChromosomes = NULL,
  pairedEnd = NULL,
  fragmentLength = NULL,
  maxFragmentLength = NULL,
  minMapq = NULL,
  removeDuplicates = NULL,
  nThreads = 1,
  verbose = TRUE
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object returned by
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- bamFiles:

  Character vector with the paths of the BAM files, in the same order as
  the samples of `counts`. Default: `NULL`, the files recorded by
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
  are reused.

- binSize:

  Numeric value with the width of the bins, in base pairs. Default:
  `10000`.

- excludeRegions:

  Logical value indicating whether the bins overlapping the regions of
  `counts` must be discarded. Default: `TRUE`.

- minCount:

  Numeric value with the minimum total count required to keep a bin.
  Default: `1`.

- excludeChromosomes:

  Character vector with the chromosomes left out of the bins and of
  their library sizes, written in either naming style, as in
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).
  `character(0)` excludes nothing. Default: `NULL`, the value used at
  the counting step.

- pairedEnd:

  Logical value, or one logical value per BAM file, indicating whether
  the reads must be counted as proper pairs. Default: `NULL`, the
  layouts resolved at the counting step.

- fragmentLength:

  Numeric value with the length to which single-end reads are extended.
  Default: `NULL`, the value used at the counting step.

- maxFragmentLength:

  Numeric value with the maximum insert size accepted for a pair.
  Default: `NULL`, the value used at the counting step.

- minMapq:

  Numeric value with the minimum mapping quality of a read. Default:
  `NULL`, the value used at the counting step.

- removeDuplicates:

  Logical value indicating whether the duplicated reads must be
  discarded. Default: `NULL`, the value used at the counting step.

- nThreads:

  Number of threads. The files are cut into pieces of at most 50 Mb,
  shared among the threads. Default: `1`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

The input `RegionSetDE.counts` object with the bin counts stored as a
`RangedSummarizedExperiment` in `metadata(counts)$background`. Its
`totals` column holds the library sizes of the bins.

## Details

Bins of ten kilobases or more are wide enough that most of them carry
background reads only, and their counts therefore track the amount of
sequencing spent outside the regions of interest. Reusing the read
parameters of
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
matters here: bins counted with a different mapping quality or duplicate
policy would return factors that do not apply to the region counts. The
parameters are taken from the object unless they are given explicitly.

The bins cover every chromosome of the BAM files except those in
`excludeChromosomes`, starting at the first base and with the last bin
of each chromosome stopping at its end. Each fragment is counted once,
in the bin holding its centre, or the 5' end of the read for single-end
data, so a fragment lying across two bins is not counted twice. The
whole genome is read whatever the value of `fullLibrarySize` used for
the regions.

## See also

[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
# The example counts already carry their background bins
counts <- loadExampleData("counts", verbose = FALSE)

backgroundBins <- S4Vectors::metadata(counts)$background
backgroundBins
#> class: RangedSummarizedExperiment 
#> dim: 1579 4 
#> metadata(4): spacing width shift bin
#> assays(1): counts
#> rownames: NULL
#> rowData names(0):
#> colnames(4): lv-H3K4me3-BN-female-bio1-tech1
#>   lv-H3K4me3-BN-male-bio2-tech1 lv-H3K4me3-SHR-male-bio2-tech1
#>   lv-H3K4me3-SHR-male-bio3-tech1
#> colData names(2): bam.files totals

# Counting them needs the BAM files the object was counted from
sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)

arCounts <- countReads(peakRegions, sampleSheet = sampleSheet, verbose = FALSE)
arCounts <- countBackground(arCounts, binSize = 10000, verbose = FALSE)

dim(S4Vectors::metadata(arCounts)$background)
#> [1] 1102    9
```
