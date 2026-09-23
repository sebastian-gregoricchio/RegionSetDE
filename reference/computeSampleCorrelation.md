# computeSampleCorrelation

Computes the correlation between the samples of an object on the log2
signal of its regions, normalised or raw, and returns it together with
the annotation of the samples, ready for
[`plotSampleCorrelation`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSampleCorrelation.md)
or for any other use.

## Usage

``` r
computeSampleCorrelation(
  object,
  set = NULL,
  contrast = NULL,
  method = "spearman",
  useOffsets = TRUE,
  topRegions = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  `RegionSetDE.counts`, `RegionSetDE.fit` or any result object of the
  package.

- set:

  Character vector with the names of the region sets used. Default:
  `NULL`, all of them.

- contrast:

  String with the name of a contrast, or its position, when `object`
  holds several of them. Default: `NULL`.

- method:

  String with the correlation, one of `"spearman"`, `"pearson"` and
  `"kendall"`. Default: `"spearman"`.

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
samples; `samples`, the `colData` of the object with a `sample` column
in the order of the matrix; and `parameters`, with the method, the
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
```
