# scoreTable

Returns the per-library table behind a `RegionSetDE.setScores` object,
one row per region set per library, with the summarised signal, the
reference it was divided by and the score.

## Usage

``` r
scoreTable(object)

# S4 method for class 'RegionSetDE.setScores'
scoreTable(object)
```

## Arguments

- object:

  `RegionSetDE.setScores` object returned by
  [`scoreRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreRegionSets.md).

## Value

A data.frame with one row per library per region set.

## See also

[`scoreRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreRegionSets.md),
[`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
setScores <- scoreRegionSets(counts, verbose = FALSE)

head(scoreTable(setScores))
#>                            sample region.set n.regions  set.signal
#> 1 lv-H3K4me3-BN-female-bio1-tech1   geneBody      1370 0.006523358
#> 2   lv-H3K4me3-BN-male-bio2-tech1   geneBody      1370 0.007104380
#> 3  lv-H3K4me3-SHR-male-bio2-tech1   geneBody      1370 0.006563504
#> 4  lv-H3K4me3-SHR-male-bio3-tech1   geneBody      1370 0.019312409
#> 5 lv-H3K4me3-BN-female-bio1-tech1 intergenic      1112 0.004901978
#> 6   lv-H3K4me3-BN-male-bio2-tech1 intergenic      1112 0.005135791
#>   reference.signal       score
#> 1      0.005052627  0.36858114
#> 2      0.005285830  0.42657861
#> 3      0.004592457  0.51519989
#> 4      0.012350258  0.64498698
#> 5      0.005052627 -0.04366946
#> 6      0.005285830 -0.04154346

# The comparisons between the sets come out of resultsTable() instead
resultsTable(setScores)
#>         set.1          set.2 n.libraries mean.delta.score   CI.lower   CI.upper
#> 1 promoterCpG promoterNonCpG           4        3.1794112  2.8832330  3.4755895
#> 2  intergenic    promoterCpG           4       -5.8874452 -6.7375513 -5.0373390
#> 3    geneBody    promoterCpG           4       -5.3558157 -6.3772026 -4.3344288
#> 4  intergenic promoterNonCpG           4       -2.7080339 -3.2874330 -2.1286349
#> 5    geneBody promoterNonCpG           4       -2.1764045 -2.9185224 -1.4342865
#> 6    geneBody     intergenic           4        0.5316295  0.3477381  0.7155209
#>   t.statistic df      p.value          FDR
#> 1   34.162894  3 5.514026e-05 0.0003308416
#> 2  -22.040162  3 2.044642e-04 0.0006133925
#> 3  -16.687697  3 4.684848e-04 0.0009369697
#> 4  -14.874329  3 6.593815e-04 0.0009890723
#> 5   -9.333139  3 2.604497e-03 0.0027156657
#> 6    9.200442  3 2.715666e-03 0.0027156657
```
