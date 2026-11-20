# Contrast accessors

Return how a result was obtained: the coefficients of the contrast, the
two groups it compares, the cut-offs its labels were set at, and, for
tiled counts, how the tiles were combined into their region.

## Usage

``` r
contrastVector(results)

# S4 method for class 'RegionSetDE.results'
contrastVector(results)

contrastGroups(results)

# S4 method for class 'RegionSetDE.results'
contrastGroups(results)

# S4 method for class 'RegionSetDE.setResults'
contrastGroups(results)

thresholdInfo(results)

# S4 method for class 'RegionSetDE.results'
thresholdInfo(results)

# S4 method for class 'RegionSetDE.setResults'
thresholdInfo(results)

# S4 method for class 'RegionSetDE.geneResults'
thresholdInfo(results)

tileCombination(results)

# S4 method for class 'RegionSetDE.results'
tileCombination(results)
```

## Arguments

- results:

  `RegionSetDE.results` object returned by
  [`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md).
  `contrastGroups` also takes a `RegionSetDE.setResults` object, and
  `thresholdInfo` a `RegionSetDE.setResults` or
  `RegionSetDE.geneResults` object.

## Value

`contrastVector` returns the coefficients of the contrast over the
columns of the design. `contrastGroups` returns a list with the `column`
of the `colData` the contrast separates and the two `groups` it
compares, the first being the level the fold change is positive for, and
for a per-region result `n.samples`, the number of samples in each; the
list is empty when the contrast is not a difference between two levels
of one variable. `thresholdInfo` returns a named list with the cut-offs,
such as `FDR` and `log2FC`, and the correction method, `adjust.method`.
`tileCombination` returns a list with the `method` that combined the
tiles of each region, empty when the counts were not tiled.

## See also

[`contrastName`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/contrastName.md),
[`contrastInfo`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/contrastInfo.md),
[`updateThresholds`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/updateThresholds.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

contrastVector(results)
#>  (Intercept) conditionSHR 
#>            0            1 
contrastGroups(results)
#> $column
#> [1] "condition"
#> 
#> $groups
#> [1] "SHR" "BN" 
#> 
#> $n.samples
#> SHR  BN 
#>   2   2 
#> 
thresholdInfo(results)
#> $FDR
#> [1] 0.05
#> 
#> $log2FC
#> [1] 0
#> 
#> $lfcThreshold
#> [1] 0
#> 
#> $adjust.method
#> [1] "BH"
#> 
tileCombination(results)
#> $applied
#> [1] FALSE
#> 
#> $method
#> [1] NA
#> 
```
