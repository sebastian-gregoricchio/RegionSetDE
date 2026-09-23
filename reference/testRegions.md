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
  signalBy = NULL,
  extraColumns = TRUE,
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

  String with the method used to combine the tiles into their region.
  `"simes"`, through
  [`csaw::combineTests`](https://rdrr.io/pkg/csaw/man/combineTests.html),
  asks whether any part of the region changes, and a single strong tile
  is enough. `"holm-min"`, through
  [`csaw::minimalTests`](https://rdrr.io/pkg/csaw/man/minimalTests.html),
  asks for several tiles to change together, three of them or 40% of the
  region when that is more, and all of them in a region shorter than
  three tiles, which suits broad domains where one tile moving on its
  own is more likely noise than biology. Default: `"simes"`.

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

- signalBy:

  String with a column of the `colData`: every level of it gets an
  `average.signal.<level>` column in the results, the average signal
  over the samples of that level alone. `FALSE` adds none. Default:
  `NULL`, the column the contrast compares two levels of, which is known
  for the three-element form and for most coefficients and expressions,
  and none when it is not.

- extraColumns:

  Annotation carried by the regions that must be appended to the result,
  at the end of the table. Either `TRUE` for every column of the
  `rowData` beyond the ones the package writes itself, `FALSE` for none,
  or a character vector naming the ones wanted. Default: `TRUE`.

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

A `RegionSetDE.results` object. Its table, read with
[`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md),
holds one row per region with the coordinates, then:

- `log2FC`: the log2 fold change, first group of the contrast over the
  second.

- `average.signal`: the average abundance of the region over all the
  samples, on the scale of the engine: the log2 counts per million of
  [`edgeR::aveLogCPM`](https://rdrr.io/pkg/edgeR/man/aveLogCPM.html) for
  edgeR, the average of the log2 values the linear model was fitted on
  for voom, limma and dream, `log2(baseMean + 1)` for DESeq2.

- `average.signal.<level>`: the same quantity computed on the samples of
  one level of `signalBy` only, one column per level.

- `stat`: the test statistic of the engine: the quasi-likelihood F of
  edgeR, or its likelihood ratio when the dispersion was fixed, the
  moderated t of voom, limma and dream, the Wald statistic of DESeq2.

- `stat.distribution`: the distribution `stat` follows under the null
  hypothesis, one among `"f"`, `"chisq"`, `"t"` and `"norm"`. `NA` for
  the threshold tests run when `lfcThreshold > 0`, whose null is not
  centred on zero.

- `df1`, `df2`: the degrees of freedom of that distribution. For the F
  of edgeR, `df1` is the numerator, the number of coefficients tested,
  and `df2` the denominator, the residual degrees of freedom plus the
  prior ones. For the likelihood ratio `df1` is the number of
  coefficients tested. For the moderated t `df1` is the total degrees of
  freedom, residual plus prior. The normal distribution of DESeq2 has
  none. `NA` where the distribution does not use them.

- `p.value`, `FDR`: the p-value and its adjustment over all the rows of
  the contrast.

- `diff.status`: `"up"`, `"down"` or `"null"`, from `FDR` and `log2FC`.

On a tiled object the statistics, the degrees of freedom and the
averages come from the tile carrying the p-value of the region, followed
by the columns the combination adds.

## Details

The multiple testing correction is applied over all the rows of the
object, across the region sets, and `regionSets` subsets the output
afterwards. Correcting inside each set separately would make the FDR of
a set depend on how many other sets were loaded, which is not a property
anyone wants in a result.

Two things follow from the combination step. With the default
`combineMethod` the p-value of a tiled region is a Simes combination, so
it answers "does any part of this region change" rather than "does the
whole region change", and a long domain that moves over one tile out of
forty will come out with a small p-value and a small overall fold
change. The `log2FC` reported for a combined region is the fold change
of the most significant tile, not an average, which is the quantity that
matches the p-value. The tile level table stays available in the `tiles`
slot, and
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

Whatever the regions were loaded with travels through to the result. A
gene name, a peak score or any other column attached to the `rowData`
comes out at the end of the table, which is what makes
[`topRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/topRegions.md)
readable and lets `plotVolcano(labelColumn = )` label the points with
something other than an identifier. On a tiled object the value is read
off the tile the combination reported, the same one the fold change
comes from, so a row describes one place rather than an average over
several.

The statistic, the distribution it follows and its degrees of freedom
are in the table so that the test can be taken further, into a power or
sample size analysis for instance.
[`contrastInfo`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/contrastInfo.md)
gathers the rest of what such an analysis needs, the engine and the
number of samples in each group of the contrast, which is also stored in
the `n.samples` element of the `contrast.groups` slot.

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

# A coefficient of the design
res <- testRegions(fit, contrast = "conditionSHR", verbose = FALSE)

# Two levels of a column, whichever of them the design took as reference
res <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
res
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

# The opposite direction, the fold changes change sign
resReverse <- testRegions(fit, contrast = c("condition", "BN", "SHR"), verbose = FALSE)

# Several contrasts on the same fit
resList <- testRegions(fit, contrast = list(shr = c("condition", "SHR", "BN"),
                                            bn = c("condition", "BN", "SHR")),
                       verbose = FALSE)
resList
#> An object of class 'RegionSetDE.resultsList'
#>   contrasts       : 2 
#> 
#>  name             contrast n.regions up down
#>   shr condition: SHR vs BN      1895  3    9
#>    bn condition: BN vs SHR      1895  9    3
topRegions(resList, contrast = "shr", FDR = 1, n = 3)
#>       region.set    region.id tile.id seqnames    start      end width
#> 1 promoterNonCpG region_02996      NA    chr12 36842295 36843294  1000
#> 2     intergenic region_03590      NA    chr12 44174500 44175499  1000
#> 3 promoterNonCpG region_00212      NA    chr12  2500829  2501828  1000
#>      log2FC average.signal average.signal.BN average.signal.SHR     stat
#> 1 -5.203835       5.159079          6.078432           2.651961 84.13617
#> 2 -2.887100       5.237329          5.930447           4.057310 43.04577
#> 3 -3.222630       4.816977          5.557234           3.451580 34.29977
#>   stat.distribution df1      df2      p.value          FDR diff.status
#> 1                 f   1 20.32506 1.148685e-08 2.176757e-05        down
#> 2                 f   1 19.05860 2.743053e-06 2.599042e-03        down
#> 3                 f   1 18.45056 1.370844e-05 6.573186e-03        down
#>       regionId
#> 1 region_02996
#> 2 region_03590
#> 3 region_00212

# Threshold inside the test rather than on the output
resStrict <- testRegions(fit, contrast = c("condition", "SHR", "BN"), lfcThreshold = 1, verbose = FALSE)
```
