# computeSampleCorrelation

Computes the correlation between the samples of an object on the log2
signal of its regions, normalised or raw, and returns it together with
the annotation of the samples, ready for
[`plotSampleCorrelation`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSampleCorrelation.md)
or for any other use. With `method = "jaccard"` the samples are compared
on where their peaks were called rather than on their signal, as the
Jaccard index of their peak sets.

## Usage

``` r
computeSampleCorrelation(
  object,
  set = NULL,
  contrast = NULL,
  samples = NULL,
  method = "spearman",
  jaccardLevel = "region",
  useOffsets = TRUE,
  topRegions = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  `RegionSetDE.counts`, `RegionSetDE.fit` or any result object of the
  package. For `method = "jaccard"`, the `RegionSetDE` object returned
  by
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md),
  or a named `GRangesList` or list of `GRanges` with the peaks of every
  sample.

- set:

  Character vector with the names of the region sets used. For
  `method = "jaccard"`, only the peaks overlapping the regions of these
  sets are compared. Default: `NULL`, all of them.

- contrast:

  String with the name of a contrast, or its position, when `object`
  holds several of them. Default: `NULL`.

- samples:

  Samples the correlation is computed on: a character vector with their
  names, a numeric vector with their positions, or a logical vector with
  one value per sample. The normalisation stored in the object is kept,
  so the values are the ones of the whole analysis restricted to these
  samples. For `method = "jaccard"` the samples whose peaks are
  compared. Default: `NULL`, all of them.

- method:

  String with the measure, one of `"spearman"`, `"pearson"` and
  `"kendall"`, computed on the signal, or `"jaccard"`, computed on the
  peak calls. Default: `"spearman"`.

- jaccardLevel:

  String with what the Jaccard index counts, either `"region"`, the
  regions of the union of all the peaks occupied by both samples over
  those occupied by either, or `"basepair"`, the base pairs covered by
  the peaks of both samples over those covered by either. Only for
  `method = "jaccard"`. Default: `"region"`.

- useOffsets:

  Logical value to indicate whether the normalisation stored in the
  object must be applied, `FALSE` scaling the samples by their library
  sizes alone. Default: `TRUE`.

- topRegions:

  Numeric value with the number of most variable regions the correlation
  is computed on. Default: `NULL`, all of them.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A list with `correlation`, the matrix of the correlations between the
samples, or of the Jaccard indices; `samples`, the `colData` of the
object with a `sample` column in the order of the matrix, or the sample
sheet the peaks were read from; and `parameters`, with the method, the
normalisation, the region sets and the number of regions used.

## Details

The values are log2 counts per million, with a prior count of 2 added to
every count and, when `useOffsets = TRUE`, the scaling factors or the
offsets stored by
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)
applied. An object that was never normalised is scaled by the library
sizes whatever `useOffsets` says, and a message reports it. The most
variable regions are chosen on the normalised values, so the raw and the
normalised correlations of one object are computed on the same rows.

A correlation does not see a single factor per sample. Scaling a library
moves its log values by a constant, which leaves the three correlations
exactly where they were, so `useOffsets` changes the matrix only when
the normalisation holds one offset per region, as `method = "loess"` and
offsets supplied from outside do. An ordination is a different matter,
and
[`plotRegionPCA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegionPCA.md)
does move with the factors.

The Jaccard index is the occupancy counterpart of the correlation, what
DiffBind draws before counting any read. It goes from zero, for two
samples with no peak in common, to one, for two identical peak sets, and
one minus it is the Jaccard distance the samples are clustered on by
[`plotSampleCorrelation`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSampleCorrelation.md).
At the `"region"` level the union of the peaks of all the samples is the
list of places, and each sample either has a peak on a place or not, so
that a broad peak and a narrow one on the same site count as a match. At
the `"basepair"` level the widths count, as in `bedtools jaccard`, and
two samples calling the same sites with different widths come out less
similar. The peaks are those kept by
[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md)
after its `blacklist` and `greylist`, so the artefacts weigh on neither
sample. The index says nothing about how strong the signal is, and a
sample with far fewer peaks than the others is dissimilar to all of them
even when its peaks are a subset of theirs: the number of peaks of every
sample, in `consensusData(x)$samples`, is worth reading next to it.
`useOffsets` and `topRegions` do not apply.

## See also

[`plotSampleCorrelation`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSampleCorrelation.md),
[`computeSamplePCA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeSamplePCA.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes

sampleCorrelation <- computeSampleCorrelation(counts, method = "pearson")
round(sampleCorrelation$correlation, 3)
#>                                 lv-H3K4me3-BN-female-bio1-tech1
#> lv-H3K4me3-BN-female-bio1-tech1                           1.000
#> lv-H3K4me3-BN-male-bio2-tech1                             0.966
#> lv-H3K4me3-SHR-male-bio2-tech1                            0.959
#> lv-H3K4me3-SHR-male-bio3-tech1                            0.926
#>                                 lv-H3K4me3-BN-male-bio2-tech1
#> lv-H3K4me3-BN-female-bio1-tech1                         0.966
#> lv-H3K4me3-BN-male-bio2-tech1                           1.000
#> lv-H3K4me3-SHR-male-bio2-tech1                          0.958
#> lv-H3K4me3-SHR-male-bio3-tech1                          0.931
#>                                 lv-H3K4me3-SHR-male-bio2-tech1
#> lv-H3K4me3-BN-female-bio1-tech1                          0.959
#> lv-H3K4me3-BN-male-bio2-tech1                            0.958
#> lv-H3K4me3-SHR-male-bio2-tech1                           1.000
#> lv-H3K4me3-SHR-male-bio3-tech1                           0.937
#>                                 lv-H3K4me3-SHR-male-bio3-tech1
#> lv-H3K4me3-BN-female-bio1-tech1                          0.926
#> lv-H3K4me3-BN-male-bio2-tech1                            0.931
#> lv-H3K4me3-SHR-male-bio2-tech1                           0.937
#> lv-H3K4me3-SHR-male-bio3-tech1                           1.000

# The samples of one strain only
brownNorway <- SummarizedExperiment::colData(counts)$condition == "BN"
round(computeSampleCorrelation(counts, samples = brownNorway, method = "pearson")$correlation, 3)
#>                                 lv-H3K4me3-BN-female-bio1-tech1
#> lv-H3K4me3-BN-female-bio1-tech1                           1.000
#> lv-H3K4me3-BN-male-bio2-tech1                             0.966
#>                                 lv-H3K4me3-BN-male-bio2-tech1
#> lv-H3K4me3-BN-female-bio1-tech1                         0.966
#> lv-H3K4me3-BN-male-bio2-tech1                           1.000

# Occupancy: the Jaccard index of the peak sets of the AR example
if (requireNamespace("consensusRegions", quietly = TRUE)) {
  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
  consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)

  peakJaccard <- computeSampleCorrelation(consensus, method = "jaccard")
  round(peakJaccard$correlation, 2)
}
#>                 AR_DMSO_r1 AR_DMSO_r2 AR_DMSO_r3 AR_R1881_4h_r1 AR_R1881_4h_r2
#> AR_DMSO_r1            1.00       0.62       0.53           0.18           0.15
#> AR_DMSO_r2            0.62       1.00       0.73           0.19           0.15
#> AR_DMSO_r3            0.53       0.73       1.00           0.18           0.13
#> AR_R1881_4h_r1        0.18       0.19       0.18           1.00           0.71
#> AR_R1881_4h_r2        0.15       0.15       0.13           0.71           1.00
#> AR_R1881_4h_r3        0.17       0.17       0.16           0.69           0.68
#> AR_R1881_24h_r1       0.12       0.12       0.12           0.54           0.65
#> AR_R1881_24h_r2       0.10       0.10       0.10           0.50           0.59
#> AR_R1881_24h_r3       0.12       0.12       0.12           0.50           0.59
#>                 AR_R1881_4h_r3 AR_R1881_24h_r1 AR_R1881_24h_r2 AR_R1881_24h_r3
#> AR_DMSO_r1                0.17            0.12            0.10            0.12
#> AR_DMSO_r2                0.17            0.12            0.10            0.12
#> AR_DMSO_r3                0.16            0.12            0.10            0.12
#> AR_R1881_4h_r1            0.69            0.54            0.50            0.50
#> AR_R1881_4h_r2            0.68            0.65            0.59            0.59
#> AR_R1881_4h_r3            1.00            0.55            0.50            0.54
#> AR_R1881_24h_r1           0.55            1.00            0.72            0.66
#> AR_R1881_24h_r2           0.50            0.72            1.00            0.66
#> AR_R1881_24h_r3           0.54            0.66            0.66            1.00
```
