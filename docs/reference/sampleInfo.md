# sampleInfo

Returns the sample table of an object as a data.frame: one row per
sample, with the annotation of the sample sheet and what the package
added along the way, such as the library sizes, the layout of each
library and, once the counts are normalised, the scaling factors. It is
[`SummarizedExperiment::colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
in a form `dplyr` and `ggplot2` take directly.

## Usage

``` r
sampleInfo(object, columns = NULL)
```

## Arguments

- object:

  `RegionSetDE.counts` or `RegionSetDE.fit` object, or any of the
  results classes, whose carried counts are used.

- columns:

  Character vector with the columns to keep, in that order. Default:
  `NULL`, all of them.

## Value

A data.frame with one row per sample.

## See also

[`libInfo`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/libInfo.md),
[`contrastInfo`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/contrastInfo.md),
[`countTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countTable.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
sampleInfo(counts)
#>                            sample
#> 1 lv-H3K4me3-BN-female-bio1-tech1
#> 2   lv-H3K4me3-BN-male-bio2-tech1
#> 3  lv-H3K4me3-SHR-male-bio2-tech1
#> 4  lv-H3K4me3-SHR-male-bio3-tech1
#>                                                                                                                     bam.file
#> 1 /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-BN-female-bio1-tech1.bam
#> 2   /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-BN-male-bio2-tech1.bam
#> 3  /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-SHR-male-bio2-tech1.bam
#> 4  /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-SHR-male-bio3-tech1.bam
#>   condition    sex biologicalReplicate paired.end library.size
#> 1        BN female                bio1      FALSE       386378
#> 2        BN   male                bio2      FALSE       400384
#> 3       SHR   male                bio2      FALSE       337600
#> 4       SHR   male                bio3      FALSE       780222

counts <- normalizeCounts(counts, method = "TMM", verbose = FALSE)
sampleInfo(counts, columns = c("sample", "condition", "library.size", "scaling.factor"))
#>                            sample condition library.size scaling.factor
#> 1 lv-H3K4me3-BN-female-bio1-tech1        BN       386378      0.7502489
#> 2   lv-H3K4me3-BN-male-bio2-tech1        BN       400384      0.7881769
#> 3  lv-H3K4me3-SHR-male-bio2-tech1       SHR       337600      0.6487672
#> 4  lv-H3K4me3-SHR-male-bio3-tech1       SHR       780222      1.8128070
```
