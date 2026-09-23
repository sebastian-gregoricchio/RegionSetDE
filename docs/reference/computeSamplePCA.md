# computeSamplePCA

Computes the principal components of the samples of an object on the
log2 signal of its most variable regions, normalised or raw, and returns
the coordinates, the variance explained and the loadings, ready for
[`plotRegionPCA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegionPCA.md)
or for any other use.

## Usage

``` r
computeSamplePCA(
  object,
  set = NULL,
  contrast = NULL,
  useOffsets = TRUE,
  topRegions = 2000,
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

- useOffsets:

  Logical value to indicate whether the normalisation stored in the
  object must be applied, `FALSE` scaling the samples by their library
  sizes alone. Default: `TRUE`.

- topRegions:

  Numeric value with the number of most variable regions the ordination
  is computed on. Default: `2000`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A list with `scores`, a data.frame with one row per sample, its
coordinates on every component and its `colData`; `variance`, a
data.frame with the standard deviation, the percentage of variance and
the cumulative percentage of every component; `loadings`, the matrix of
the weights of the regions on the components; and `parameters`, with the
region sets, the normalisation and the number of regions used.

## Details

The values are log2 counts per million, with a prior count of 2 added to
every count and, when `useOffsets = TRUE`, the scaling factors or the
offsets stored by
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)
applied. The regions are centred and not scaled, so a region weighs as
much as it varies. An object that was never normalised is scaled by the
library sizes whatever `useOffsets` says, and a message reports it. The
most variable regions are chosen on the normalised values, so the raw
and the normalised ordinations of one object are computed on the same
rows.

## See also

[`plotRegionPCA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegionPCA.md),
[`computeSampleCorrelation`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeSampleCorrelation.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes

samplePCA <- computeSamplePCA(counts, topRegions = 1000)
samplePCA$variance
#>   component standard.deviation variance.percent cumulative.percent
#> 1       PC1          28.874387        84.205359           84.20536
#> 2       PC2           9.753035         9.607131           93.81249
#> 3       PC3           7.827100         6.187510          100.00000
head(samplePCA$scores[, c("sample", "PC1", "PC2", "condition")])
#>                            sample       PC1         PC2 condition
#> 1 lv-H3K4me3-BN-female-bio1-tech1 -17.38023   4.6003874        BN
#> 2   lv-H3K4me3-BN-male-bio2-tech1 -14.40160   8.3444835        BN
#> 3  lv-H3K4me3-SHR-male-bio2-tech1 -11.37325 -13.9150911       SHR
#> 4  lv-H3K4me3-SHR-male-bio3-tech1  43.15508   0.9702203       SHR
```
