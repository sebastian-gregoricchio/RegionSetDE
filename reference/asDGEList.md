# asDGEList

Turns a `RegionSetDE.counts` object into a `DGEList`, with the
normalisation stored in the object carried across as offsets and the
region annotation kept in the `genes` slot.

## Usage

``` r
asDGEList(counts, assay = "counts", useOffsets = TRUE, verbose = TRUE)
```

## Arguments

- counts:

  `RegionSetDE.counts` object, or a `RegionSetDE.fit`, in which case the
  list the model was fitted on is returned as it is.

- assay:

  String with the name of the assay holding the counts. Default:
  `"counts"`.

- useOffsets:

  Logical value to indicate whether the normalisation stored in the
  object must be carried across. Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `DGEList`.

## Details

This exists because the conversion is easy to get backwards. The object
stores a `scaling.factor` that normalised values are divided by, while a
generalised linear model wants an offset on the log scale of the
effective library size, and the two differ by a sign and a constant.
Writing the factors into `lib.size` by hand, or exponentiating them the
wrong way, produces a `DGEList` that fits without complaint and reports
fold changes in the wrong direction. Here the factors go through
[`edgeR::scaleOffset`](https://rdrr.io/pkg/edgeR/man/scaleOffset.html),
which puts them on the scale of the library sizes and leaves the
coefficients readable, and there is no direction left to get wrong.

`as(counts, "DGEList")` does the same with the defaults.

## See also

[`asDESeqDataSet`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/asDESeqDataSet.md),
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md),
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes

dgeList <- asDGEList(counts)
#> The normalisation has been carried across as offsets, so do not set 'norm.factors' or 'lib.size' on top of it.
dgeList$samples
#>                                 group lib.size norm.factors
#> lv-H3K4me3-BN-female-bio1-tech1     1   386378            1
#> lv-H3K4me3-BN-male-bio2-tech1       1   400384            1
#> lv-H3K4me3-SHR-male-bio2-tech1      1   337600            1
#> lv-H3K4me3-SHR-male-bio3-tech1      1   780222            1
#>                                                          sample
#> lv-H3K4me3-BN-female-bio1-tech1 lv-H3K4me3-BN-female-bio1-tech1
#> lv-H3K4me3-BN-male-bio2-tech1     lv-H3K4me3-BN-male-bio2-tech1
#> lv-H3K4me3-SHR-male-bio2-tech1   lv-H3K4me3-SHR-male-bio2-tech1
#> lv-H3K4me3-SHR-male-bio3-tech1   lv-H3K4me3-SHR-male-bio3-tech1
#>                                                                                                                                                   bam.file
#> lv-H3K4me3-BN-female-bio1-tech1 /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.4/chromstaRData/extdata/euratrans/lv-H3K4me3-BN-female-bio1-tech1.bam
#> lv-H3K4me3-BN-male-bio2-tech1     /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.4/chromstaRData/extdata/euratrans/lv-H3K4me3-BN-male-bio2-tech1.bam
#> lv-H3K4me3-SHR-male-bio2-tech1   /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.4/chromstaRData/extdata/euratrans/lv-H3K4me3-SHR-male-bio2-tech1.bam
#> lv-H3K4me3-SHR-male-bio3-tech1   /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.4/chromstaRData/extdata/euratrans/lv-H3K4me3-SHR-male-bio3-tech1.bam
#>                                 condition    sex biologicalReplicate paired.end
#> lv-H3K4me3-BN-female-bio1-tech1        BN female                bio1      FALSE
#> lv-H3K4me3-BN-male-bio2-tech1          BN   male                bio2      FALSE
#> lv-H3K4me3-SHR-male-bio2-tech1        SHR   male                bio2      FALSE
#> lv-H3K4me3-SHR-male-bio3-tech1        SHR   male                bio3      FALSE
#>                                 library.size norm.factor scaling.factor
#> lv-H3K4me3-BN-female-bio1-tech1       386378   0.7959318      0.5433126
#> lv-H3K4me3-BN-male-bio2-tech1         400384   0.8036577      0.5684724
#> lv-H3K4me3-SHR-male-bio2-tech1        337600   0.9214493      0.5495858
#> lv-H3K4me3-SHR-male-bio3-tech1        780222   1.6966084      2.3386292
```
