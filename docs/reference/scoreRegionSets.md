# scoreRegionSets

Computes one signal score per region set per library, without a
contrast, so that the sets can be compared to each other in an
experiment holding a single condition. The score is the summarised
signal over the regions of a set divided by a reference measured in the
same library, and the sets are then compared library by library, which
puts the replication in the libraries rather than in the regions. It
answers whether one set carries more signal than another. It does not
answer whether that difference comes from the factor or from what the
sets are made of.

## Usage

``` r
scoreRegionSets(
  counts,
  regionSets = NULL,
  comparisons = NULL,
  reference = "background",
  assay = NULL,
  perBasepair = TRUE,
  summary = "mean",
  pseudoCount = 0.5,
  minRegions = 10,
  confLevel = 0.95,
  adjustMethod = "BH",
  verbose = TRUE
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object returned by
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- regionSets:

  Character vector with the names of the sets to score. Default: `NULL`,
  all of them.

- comparisons:

  List of character vectors of length two, naming the sets compared to
  each other. Default: `NULL`, every pair.

- reference:

  String with what the signal of a set is divided by, one among
  `"background"` (the bins of
  [`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md)),
  `"regions"` (every region loaded in the object) and `"none"` (no
  division). Default: `"background"`.

- assay:

  String with the name of the assay holding the signal, ignored when
  `reference` is `"background"`. Default: `NULL`, the normalised assay
  when present and the raw counts otherwise.

- perBasepair:

  Logical value indicating whether the signal must be divided by the
  width of the region before it is summarised. Default: `TRUE`.

- summary:

  String with how the regions of a set are summarised, either `"mean"`
  or `"median"`. Default: `"mean"`.

- pseudoCount:

  Numeric value added to the signal before the ratio is taken. Default:
  `0.5`.

- minRegions:

  Numeric value with the number of regions a set needs to be scored.
  Default: `10`.

- confLevel:

  Numeric value with the level of the confidence interval on the
  difference between two sets. Default: `0.95`.

- adjustMethod:

  String with the multiple testing correction applied across the
  comparisons. Default: `"BH"`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.setScores` object.
[`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md)
returns the comparisons between the sets,
[`scoreTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreTable.md)
the per-library scores they were computed from, and
[`plotSetSignal`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetSignal.md)
draws both.

## Details

Every test in the package runs on a contrast, so a design holding one
condition and three replicates leaves
[`testSetContrast`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md)
with nothing to compare. This function covers that case, and the price
of covering it is worth stating before the output is read.

The score of a set is a ratio taken inside one library, so the
sequencing depth and any scaling factor shared by the numerator and the
denominator cancel before the score exists. The background bins carry
raw counts, so with `reference = "background"` the regions are read from
the `counts` assay as well and the normalisation never enters, which
also means that changing it cannot move the result. With
`reference = "none"` that protection is gone: the score becomes the
signal itself, the libraries sit on their own scales again, and the
normalised assay is the one to pass.

The reference cancels a second time when two sets are compared, since
both scores of a library were divided by the same number. What the
comparison tests is the log ratio between the summarised signal of the
two sets, taken library by library. It depends on neither the reference
nor the normalisation. It does depend on the libraries, and there are
usually three of them, so the interval rests on two degrees of freedom
and comes out wide. That width is what three replicates support.

The regions are not the replication here and are deliberately not
treated as it. Thirty thousand promoters collapse to one number per
library, and the variance reaching the test is the variance between
libraries. Summarising over the regions instead of testing across them
is why the p-value is not the vanishing one a per-region test over the
same data would return.

What no arrangement of these numbers separates is the factor from the
regions. Two sets differing in width, mappability, GC content or
accessibility differ in coverage in a library where nothing is bound,
and that difference repeats across replicates in the same way real
binding does. `perBasepair` takes away the crudest part of it and leaves
the rest standing. A difference between two sets is therefore a
difference in signal, and calling it a difference in factor is a
separate argument, made with the composition of the sets shown next to
the result. An input, an IgG or a spike-in turns the score into an
enrichment and removes the need to make that argument at all; where one
exists,
[`testSetContrast`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md)
on the IP against it is the better instrument.

## See also

[`testSetContrast`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md),
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md),
[`plotSetSignal`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetSignal.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)

# One score per library per set, plus every pairwise comparison
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
```
