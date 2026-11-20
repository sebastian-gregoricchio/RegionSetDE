# RegionSetDE.counts class

S4 class storing the read counts computed over a collection of region
sets. It extends `RangedSummarizedExperiment`, therefore `assay`,
`colData`, `rowRanges` and the subsetting operators behave as usual,
while the filters applied upstream remain accessible in the inherited
provenance slots.

## Slots

- `counting.level`:

  String indicating whether the rows are regions (`"region"`) or tiles
  of a region (`"tile"`). It is the same vocabulary the
  `RegionSetDE.fit` and `RegionSetDE.results` classes use, since the
  value travels from here into both of them.

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
counts
#> class: RegionSetDE.counts 
#> dim: 3224 4 
#> metadata(3): signal.type count.like background
#> assays(1): counts
#> rownames(3224): promoterNonCpG|region_00002 promoterNonCpG|region_00003
#>   ... promoterCpG|region_03797 promoterCpG|region_03798
#> rowData names(4): region.set region.id tile.id regionId
#> colnames(4): lv-H3K4me3-BN-female-bio1-tech1
#>   lv-H3K4me3-BN-male-bio2-tech1 lv-H3K4me3-SHR-male-bio2-tech1
#>   lv-H3K4me3-SHR-male-bio3-tech1
#> colData names(7): sample bam.file ... paired.end library.size
head(sampleInfo(counts))
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
```
