# pairwiseContrasts

Writes out the contrasts between the levels of a column of the sample
table, ready to be handed to
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
or
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md):
every pair of levels, or every level against a reference.

## Usage

``` r
pairwiseContrasts(object, column, levels = NULL, reference = NULL)
```

## Arguments

- object:

  `RegionSetDE.counts` or `RegionSetDE.fit` object, whose `colData`
  holds the column. A data.frame with the sample annotation is accepted
  as well.

- column:

  String with the name of the column.

- levels:

  Character vector with the levels to use, in the order that sets the
  direction of the contrasts: each level is compared against the ones
  before it. Default: `NULL`, the levels of the column when it is a
  factor, and the order in which they first appear otherwise.

- reference:

  String with a level to compare every other level against. When given,
  only those contrasts are written, instead of every pair. Default:
  `NULL`.

## Value

A named list of character vectors `c(column, levelA, levelB)`, one per
contrast, each read as `levelA` against `levelB`. The names are
`levelA_vs_levelB`, and they become the names of the contrasts in the
results.

## Details

With *k* levels there are *k(k - 1) / 2* pairs, which grows quickly, and
every pair is a separate family of tests. The correction for multiple
testing is applied within each contrast, never across them, so the more
contrasts are tested the more of the reported regions are false
positives overall. A reference level keeps the number at *k - 1*, and is
usually what a treatment design is asking.

## See also

[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
pairwiseContrasts(counts, column = "condition")
#> $SHR_vs_BN
#> [1] "condition" "SHR"       "BN"       
#> 

sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)

# Every pair of the three conditions
pairwiseContrasts(sampleSheet, column = "condition")
#> $R1881_4h_vs_DMSO
#> [1] "condition" "R1881_4h"  "DMSO"     
#> 
#> $R1881_24h_vs_DMSO
#> [1] "condition" "R1881_24h" "DMSO"     
#> 
#> $R1881_24h_vs_R1881_4h
#> [1] "condition" "R1881_24h" "R1881_4h" 
#> 

# Every treatment against the vehicle only
pairwiseContrasts(sampleSheet, column = "condition", reference = "DMSO")
#> $R1881_4h_vs_DMSO
#> [1] "condition" "R1881_4h"  "DMSO"     
#> 
#> $R1881_24h_vs_DMSO
#> [1] "condition" "R1881_24h" "DMSO"     
#> 

# The list goes straight into testRegions()
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = pairwiseContrasts(fit, column = "condition"), verbose = FALSE)
results
#> An object of class 'RegionSetDE.resultsList'
#>   contrasts       : 1 
#> 
#>       name             contrast n.regions up down
#>  SHR_vs_BN condition: SHR vs BN      1895  3    9
```
