# resultsTable

Returns the per-region table of a `RegionSetDE.results` object, or the
tables of every contrast of a list of them, stacked.

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

# S4 method for class 'RegionSetDE.setScores'
resultsTable(results)
```

## Arguments

- results:

  `RegionSetDE.results` object, or one of the list classes holding
  several contrasts. A single contrast is taken out of a list by name,
  `results$name` or `results[["name"]]`.

## Value

A data.frame with one row per region. For a list of contrasts the tables
are stacked, and the first two columns say which contrast a row belongs
to: `contrast`, the name the contrast was given, which is the string the
`contrast` argument of
[`topRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/topRegions.md),
[`plotVolcano`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotVolcano.md)
and the other functions takes, and `contrast.description`, what it
compares.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

resultTable <- resultsTable(results)
head(resultTable)
#>       region.set    region.id tile.id seqnames  start    end width     log2FC
#> 1 promoterNonCpG region_00012      NA    chr12  26988  27987  1000 -0.4744875
#> 2 promoterNonCpG region_00017      NA    chr12  39449  40448  1000 -1.6911510
#> 3 promoterNonCpG region_00019      NA    chr12  44116  45115  1000  0.1328365
#> 4 promoterNonCpG region_00020      NA    chr12  46527  47526  1000  1.3000608
#> 5 promoterNonCpG region_00026      NA    chr12  89229  90228  1000 -0.0791789
#> 6 promoterNonCpG region_00043      NA    chr12 330751 331750  1000 -1.2370885
#>   average.signal average.signal.BN average.signal.SHR        stat
#> 1       3.098695          3.154609           3.071096 0.166290413
#> 2       3.141973          3.512708           2.774385 2.337588301
#> 3       3.273341          3.149796           3.390540 0.010028735
#> 4       3.490690          2.930994           3.854730 2.084865356
#> 5       3.220877          3.152202           3.293212 0.007114484
#> 6       3.411317          3.668027           3.170663 2.273415136
#>   stat.distribution df1      df2   p.value       FDR diff.status     regionId
#> 1                 f   1 18.94543 0.6880013 0.9200864        null region_00012
#> 2                 f   1 18.73556 0.1429932 0.7997195        null region_00017
#> 3                 f   1 18.45274 0.9213104 0.9807286        null region_00019
#> 4                 f   1 20.64203 0.1637725 0.7997195        null region_00020
#> 5                 f   1 20.56519 0.9335959 0.9807286        null region_00026
#> 6                 f   1 20.04174 0.1472119 0.7997195        null region_00043

table(resultTable$region.set, resultTable$diff.status)
#>                 
#>                  down null  up
#>   geneBody          3  906   0
#>   intergenic        3  435   2
#>   promoterCpG       1  268   0
#>   promoterNonCpG    2  274   1

```
