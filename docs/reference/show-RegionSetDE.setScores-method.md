# show method for RegionSetDE.setScores

Prints a summary of the region set scores and of the comparisons between
them.

## Usage

``` r
# S4 method for class 'RegionSetDE.setScores'
show(object)
```

## Arguments

- object:

  `RegionSetDE.setScores` object.

## Value

Prints the summary to the console.

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
scoreRegionSets(counts, verbose = FALSE)
#> An object of class 'RegionSetDE.setScores'
#>   region sets     : geneBody, intergenic, promoterCpG, promoterNonCpG 
#>   libraries       : 4 
#>   reference       : background 
#>   signal          : mean of 'counts' per base pair 
#>   comparisons     : 6 
#> 
#>        set.1          set.2 n.libraries mean.delta.score CI.lower CI.upper
#>  promoterCpG promoterNonCpG           4            3.180    2.880    3.480
#>   intergenic    promoterCpG           4           -5.890   -6.740   -5.040
#>     geneBody    promoterCpG           4           -5.360   -6.380   -4.330
#>   intergenic promoterNonCpG           4           -2.710   -3.290   -2.130
#>     geneBody promoterNonCpG           4           -2.180   -2.920   -1.430
#>     geneBody     intergenic           4            0.532    0.348    0.716
#>   p.value      FDR
#>  5.51e-05 0.000331
#>  2.04e-04 0.000613
#>  4.68e-04 0.000937
#>  6.59e-04 0.000989
#>  2.60e-03 0.002720
#>  2.72e-03 0.002720
#> 
#> The libraries are the replication, the composition of the sets is not controlled for.
```
