# Set score accessors

Return the parts of a `RegionSetDE.setScores` object beyond the
per-library scores of
[`scoreTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreTable.md):
the comparisons between the sets, the samples behind the libraries, and
how the signal of a set was summarised.

## Usage

``` r
scoreComparisons(scores)

# S4 method for class 'RegionSetDE.setScores'
scoreComparisons(scores)

scoreSamples(scores)

# S4 method for class 'RegionSetDE.setScores'
scoreSamples(scores)

scoreReference(scores)

# S4 method for class 'RegionSetDE.setScores'
scoreReference(scores)

scoreAssay(scores)

# S4 method for class 'RegionSetDE.setScores'
scoreAssay(scores)

scoreSummary(scores)

# S4 method for class 'RegionSetDE.setScores'
scoreSummary(scores)

scorePerBasepair(scores)

# S4 method for class 'RegionSetDE.setScores'
scorePerBasepair(scores)
```

## Arguments

- scores:

  `RegionSetDE.setScores` object returned by
  [`scoreRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreRegionSets.md).

## Value

`scoreComparisons` returns a data.frame with one row per pair of sets,
as `resultsTable` does. `scoreSamples` returns a data.frame with the
`colData` of the counts the scores came from. `scoreReference` returns
what the signal of a set was divided by, one among `"background"`,
`"regions"` and `"none"`. `scoreAssay` returns the name of the assay the
signal was read from. `scoreSummary` returns `"mean"` or `"median"`, the
summary of the regions of a set. `scorePerBasepair` returns `TRUE` when
the signal was divided by the width of the regions before being
summarised.

## See also

[`scoreTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreTable.md),
[`scoreRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreRegionSets.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
setScores <- scoreRegionSets(counts, verbose = FALSE)

scoreComparisons(setScores)
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
head(scoreSamples(setScores), 3)
#>                            sample                        sample.1
#> 1 lv-H3K4me3-BN-female-bio1-tech1 lv-H3K4me3-BN-female-bio1-tech1
#> 2   lv-H3K4me3-BN-male-bio2-tech1   lv-H3K4me3-BN-male-bio2-tech1
#> 3  lv-H3K4me3-SHR-male-bio2-tech1  lv-H3K4me3-SHR-male-bio2-tech1
#>                                                                                                                     bam.file
#> 1 /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-BN-female-bio1-tech1.bam
#> 2   /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-BN-male-bio2-tech1.bam
#> 3  /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-SHR-male-bio2-tech1.bam
#>   condition    sex biologicalReplicate paired.end library.size
#> 1        BN female                bio1      FALSE       386378
#> 2        BN   male                bio2      FALSE       400384
#> 3       SHR   male                bio2      FALSE       337600
scoreReference(setScores)
#> [1] "background"
scoreAssay(setScores)
#> [1] "counts"
scoreSummary(setScores)
#> [1] "mean"
scorePerBasepair(setScores)
#> [1] TRUE
```
