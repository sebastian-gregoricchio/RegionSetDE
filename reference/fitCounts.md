# fitCounts

Returns the `RegionSetDE.counts` object on which a model was fitted.

## Usage

``` r
fitCounts(fit)

# S4 method for class 'RegionSetDE.fit'
fitCounts(fit)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

## Value

A `RegionSetDE.counts` object.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

# The counts the model was fitted on, normalisation and filtering included
fitCounts(fit)
#> class: RegionSetDE.counts 
#> dim: 1895 4 
#> metadata(3): signal.type background normalization
#> assays(2): counts norm.counts
#> rownames(1895): promoterNonCpG|region_00012 promoterNonCpG|region_00017
#>   ... promoterCpG|region_03797 promoterCpG|region_03798
#> rowData names(3): region.set region.id tile.id
#> colnames(4): lv-H3K4me3-BN-female-bio1-tech1
#>   lv-H3K4me3-BN-male-bio2-tech1 lv-H3K4me3-SHR-male-bio2-tech1
#>   lv-H3K4me3-SHR-male-bio3-tech1
#> colData names(9): sample bam.file ... norm.factor scaling.factor

# The underlying engine object, here an edgeR fit
class(fitObject(fit))
#> [1] "DGEGLM"
#> attr(,"package")
#> [1] "edgeR"

```
