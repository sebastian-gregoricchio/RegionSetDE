# RegionSetDE example data

Description of the example objects installed under `inst/extdata` and of
the dataset they come from. The objects are read with
[`loadExampleData`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadExampleData.md),
which is where the individual formats are documented.

## See also

[`loadExampleData`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadExampleData.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
# The counts of the EURATRANS example and the sheet of the AR example
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

sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
head(sampleSheet[, c("sample", "condition")])
#>           sample condition
#> 1     AR_DMSO_r1      DMSO
#> 2     AR_DMSO_r2      DMSO
#> 3     AR_DMSO_r3      DMSO
#> 4 AR_R1881_4h_r1  R1881_4h
#> 5 AR_R1881_4h_r2  R1881_4h
#> 6 AR_R1881_4h_r3  R1881_4h
```
