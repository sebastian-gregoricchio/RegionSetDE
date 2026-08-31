# testRegions

Tests a contrast on a `RegionSetDE.fit` object and returns one row per
region. When the counts were tiled, every tile is tested on its own and
the p-values are then combined back to the region, so that the region
stays the unit of inference even though the signal was measured at a
finer scale.

## Usage

``` r
testRegions(
  fit,
  contrast,
  combine = TRUE,
  combineMethod = "simes",
  lfcThreshold = 0,
  FDR = 0.05,
  log2FC = 0,
  adjustMethod = "BH",
  regionSets = NULL,
  carryCounts = TRUE,
  verbose = TRUE
)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

- contrast:

  Contrast to test, given in one of four ways. A character vector of
  length three, `c("column", "groupA", "groupB")`, naming a column of
  the `colData` and two of its levels, which is the form to reach for
  when the design uses a reference level. A string with the name of a
  design column, e.g. `"conditionCOMBO"`. A string written as an
  expression over the design columns, e.g.
  `"conditionCOMBO - conditionEPZ"`. Or a numeric vector with one
  coefficient per column of the design. A named list of any of these
  runs every contrast on the same fit and returns a
  `RegionSetDE.resultsList`.

- combine:

  Logical value to indicate whether the tile level p-values must be
  combined into one value per region. Ignored when the counts were not
  tiled. Default: `TRUE`.

- combineMethod:

  String with the method used to combine the tiles, among those accepted
  by
  [`csaw::combineTests`](https://rdrr.io/pkg/csaw/man/combineTests.html):
  `"simes"`, `"holm-min"`, `"wilcoxon"` and `"stouffer"`. Default:
  `"simes"`.

- lfcThreshold:

  Numeric value with the log2 fold change against which the null
  hypothesis is tested. A value above zero moves the threshold inside
  the test, through
  [`edgeR::glmTreat`](https://rdrr.io/pkg/edgeR/man/glmTreat.html),
  [`limma::treat`](https://rdrr.io/pkg/limma/man/ebayes.html) or the
  `lfcThreshold` of
  [`DESeq2::results`](https://rdrr.io/pkg/DESeq2/man/results.html),
  which is stricter and better calibrated than filtering the output
  afterwards. Default: `0`.

- FDR:

  Numeric value with the adjusted p-value cut-off used to fill the
  `diff.status` column. Default: `0.05`.

- log2FC:

  Numeric value with the absolute log2 fold change cut-off used to fill
  the `diff.status` column. Default: `0`.

- adjustMethod:

  String with the multiple testing correction, passed to
  [`stats::p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default:
  `"BH"`.

- regionSets:

  Character vector with the names of the region sets to keep in the
  output. Default: `NULL`, all of them.

- carryCounts:

  Logical value to indicate whether the counts must travel inside the
  result, so that
  [`plotRegion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegion.md)
  and
  [`plotTopHeatmap`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotTopHeatmap.md)
  can draw the values without being handed the counts object again.
  Several contrasts run on one fit share the same copy in memory.
  Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.results` object.

## Details

The multiple testing correction is applied over all the rows of the
object, across the region sets, and `regionSets` subsets the output
afterwards. Correcting inside each set separately would make the FDR of
a set depend on how many other sets were loaded, which is not a property
anyone wants in a result.

Two things follow from the combination step. The p-value of a tiled
region is a Simes combination, so it answers "does any part of this
region change" rather than "does the whole region change", and a long
domain that moves over one tile out of forty will come out with a small
p-value and a small overall fold change. The `log2FC` reported for a
combined region is the fold change of the most significant tile, not an
average, which is the quantity that matches the p-value. The tile level
table stays available in the `tiles` slot, and
[`plotRegion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegion.md)
draws it.

A design written as `~ condition` spends one coefficient per level
except the first, so a level can be a coefficient in the design or the
reference the others are measured against, depending on how the factor
was ordered. Naming a coefficient that turns out to be the reference is
the usual source of confusion, and it is what
`c("column", "groupA", "groupB")` avoids: that form averages the design
rows of each group and takes the difference, which gives the same
contrast whatever the reference is and whether the design was written as
`~ condition` or `~ 0 + condition`. With other covariates in the design
the averaging picks up their imbalance between the two groups, so it
describes what it says only when the design is reasonably balanced.

The `diff.status` column is a labelling convenience, not a claim. It is
filled from `FDR` and `log2FC` and used by the plotting functions; the
thresholds are stored in the object so that a figure can state them.

## See also

[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md),
[`topRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/topRegions.md),
[`plotVolcano`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotVolcano.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

# The three-element form works whatever the reference level is
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
results
#> An object of class 'RegionSetDE.results'
#>   contrast        : condition: SHR vs BN 
#>   engine          : edgeR 
#>   regions         : 1895 
#>   counts carried  : 4 samples
#>   thresholds      : FDR < 0.05 | |log2FC| > 0 
#>   changing regions:
#>     promoterNonCpG: 1 up, 2 down
#>     intergenic: 2 up, 3 down
#>     geneBody: 0 up, 3 down
#>     promoterCpG: 0 up, 1 down

head(resultsTable(results))
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

# Testing against a fold change threshold rather than against zero
strictResults <- testRegions(fit, contrast = c("condition", "SHR", "BN"),
                             lfcThreshold = 1, verbose = FALSE)
strictResults
#> An object of class 'RegionSetDE.results'
#>   contrast        : condition: SHR vs BN 
#>   engine          : edgeR 
#>   regions         : 1895 
#>   counts carried  : 4 samples
#>   thresholds      : FDR < 0.05 | |log2FC| > 0 
#>   changing regions:
#>     promoterNonCpG: 0 up, 1 down
#>     intergenic: 0 up, 1 down
#>     geneBody: 0 up, 0 down
#>     promoterCpG: 0 up, 0 down
```
