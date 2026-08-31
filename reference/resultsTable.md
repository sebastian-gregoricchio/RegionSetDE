# resultsTable

Returns the per-region table of a `RegionSetDE.results` object.

## Usage

``` r
resultsTable(results)

# S4 method for class 'RegionSetDE.results'
resultsTable(results)

# S4 method for class 'RegionSetDE.setResults'
resultsTable(results)

# S4 method for class 'RegionSetDE.resultsList'
resultsTable(results)

# S4 method for class 'RegionSetDE.setResultsList'
resultsTable(results)
```

## Arguments

- results:

  `RegionSetDE.results` object.

## Value

A data.frame with one row per region.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

resultTable <- resultsTable(results)
head(resultTable)
#>       region.set    region.id seqnames  start    end width     log2FC
#> 1 promoterNonCpG region_00012    chr12  26988  27987  1000 -0.4744875
#> 2 promoterNonCpG region_00017    chr12  39449  40448  1000 -1.6911510
#> 3 promoterNonCpG region_00019    chr12  44116  45115  1000  0.1328365
#> 4 promoterNonCpG region_00020    chr12  46527  47526  1000  1.3000608
#> 5 promoterNonCpG region_00026    chr12  89229  90228  1000 -0.0791789
#> 6 promoterNonCpG region_00043    chr12 330751 331750  1000 -1.2370885
#>   average.signal        stat   p.value       FDR diff.status
#> 1       3.098695 0.166290413 0.6880013 0.9200864        null
#> 2       3.141973 2.337588301 0.1429932 0.7997195        null
#> 3       3.273341 0.010028735 0.9213104 0.9807286        null
#> 4       3.490690 2.084865356 0.1637725 0.7997195        null
#> 5       3.220877 0.007114484 0.9335959 0.9807286        null
#> 6       3.411317 2.273415136 0.1472119 0.7997195        null

table(resultTable$region.set, resultTable$diff.status)
#>                 
#>                  down null  up
#>   geneBody          3  906   0
#>   intergenic        3  435   2
#>   promoterCpG       1  268   0
#>   promoterNonCpG    2  274   1

```
