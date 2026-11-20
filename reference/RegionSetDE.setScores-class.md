# RegionSetDE.setScores class

S4 class storing the signal scores computed by
[`scoreRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreRegionSets.md),
one per region set per library, together with the paired comparisons
between the sets. No contrast lies behind these numbers and no model was
fitted on them, which is what separates this class from
`RegionSetDE.setResults`: the replication is the libraries, and the
comparison runs between sets rather than between conditions.

## Usage

``` r
# S4 method for class 'RegionSetDE.setScores'
show(object)
```

## Arguments

- object:

  `RegionSetDE.setScores` object.

## Value

[`show()`](https://rdrr.io/r/methods/show.html) prints a summary of the
object to the console.

## Slots

- `scores`:

  `DataFrame` with one row per library per set, carrying the summarised
  signal, the reference it was divided by and the `score` on a log2
  scale.

- `comparisons`:

  `DataFrame` with one row per pair of sets, carrying the mean paired
  difference between their scores, its confidence interval, the paired
  t-statistic and the adjusted p-value.

- `sample.metadata`:

  `DataFrame` with the `colData` of the counts object the scores came
  from, so that the libraries can be grouped and coloured downstream
  without handing that object over again.

- `reference`:

  String naming what the signal of a set was divided by, one among
  `"background"`, `"regions"` and `"none"`.

- `assay`:

  String with the name of the assay the signal was read from.

- `summary`:

  String with how the regions of a set were summarised, either `"mean"`
  or `"median"`.

- `per.basepair`:

  Logical value indicating whether the signal was divided by the width
  of the region before being summarised.

## See also

[`scoreRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreRegionSets.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)

setScores <- scoreRegionSets(counts, verbose = FALSE)
setScores
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
