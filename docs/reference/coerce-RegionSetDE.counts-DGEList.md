# Coerce a counts object to a DGEList

Registers the `as(x, "DGEList")` idiom, which calls
[`asDGEList`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/asDGEList.md)
with its defaults.

## Value

A `DGEList` holding the counts, the sample metadata, the region
annotation and the normalisation as offsets.

## See also

[`asDGEList`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/asDGEList.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes

dgeList <- as(counts, "DGEList")
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
#> lv-H3K4me3-BN-female-bio1-tech1 /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-BN-female-bio1-tech1.bam
#> lv-H3K4me3-BN-male-bio2-tech1     /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-BN-male-bio2-tech1.bam
#> lv-H3K4me3-SHR-male-bio2-tech1   /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-SHR-male-bio2-tech1.bam
#> lv-H3K4me3-SHR-male-bio3-tech1   /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-SHR-male-bio3-tech1.bam
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
