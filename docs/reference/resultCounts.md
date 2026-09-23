# resultCounts

Returns the `RegionSetDE.counts` object carried inside a result, the one
the contrast was computed on.

## Usage

``` r
resultCounts(results)

# S4 method for class 'RegionSetDE.results'
resultCounts(results)

# S4 method for class 'RegionSetDE.setResults'
resultCounts(results)

# S4 method for class 'RegionSetDE.resultsList'
resultCounts(results)

# S4 method for class 'RegionSetDE.setResultsList'
resultCounts(results)
```

## Arguments

- results:

  `RegionSetDE.results` object.

## Value

A `RegionSetDE.counts` object.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

resultCounts(results)
#> class: RegionSetDE.counts 
#> dim: 1895 4 
#> metadata(5): signal.type count.like background background.holdout
#>   normalization
#> assays(2): counts norm.counts
#> rownames(1895): promoterNonCpG|region_00012 promoterNonCpG|region_00017
#>   ... promoterCpG|region_03797 promoterCpG|region_03798
#> rowData names(4): region.set region.id tile.id regionId
#> colnames(4): lv-H3K4me3-BN-female-bio1-tech1
#>   lv-H3K4me3-BN-male-bio2-tech1 lv-H3K4me3-SHR-male-bio2-tech1
#>   lv-H3K4me3-SHR-male-bio3-tech1
#> colData names(9): sample bam.file ... norm.factor scaling.factor
```
