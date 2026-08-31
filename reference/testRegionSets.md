# testRegionSets

Asks whether a region set responds to a contrast as a whole. Two
questions are answered side by side: whether the regions of the set move
away from zero, which is a self-contained claim, and whether they move
more than the regions they are compared against, which is a competitive
one. Both are computed from the per-region statistics of the same fit,
so they never disagree with
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
on the design, the offsets or the dispersion.

## Usage

``` r
testRegionSets(
  fit,
  contrast,
  method = c("camera", "fry"),
  universe = NULL,
  matchOn = c("width", "abundance"),
  universeRatio = 5,
  interRegionCor = NULL,
  useRanks = FALSE,
  FDR = 0.05,
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

  Contrast to test, in the syntax accepted by
  [`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
  or a named list of contrasts to run on the same fit.

- method:

  Character vector with the tests to run, among `"camera"` (competitive)
  and `"fry"` (self-contained). Default: `c("camera", "fry")`.

- universe:

  What each set is compared against in the competitive test. Default:
  `NULL`, the universe carried by the fit. A `RegionSetDE.universe`
  object, or the strings `"matched"` and `"all"`, override it and are
  built here.

- matchOn:

  Character vector with the covariates the comparison rows are matched
  on, when one has to be built here. Default: `c("width", "abundance")`.

- universeRatio:

  Numeric value with the number of comparison rows drawn per region of
  the set, when one has to be built here. Default: `5`.

- interRegionCor:

  Numeric value with the correlation between the regions of a set, used
  to inflate the variance. Default: `NULL`, estimated from the residuals
  of each set, or held at 0.01 when the design leaves no residual to
  estimate it from.

- useRanks:

  Logical value to indicate whether `camera` must work on the ranks
  rather than on the statistics, which is more robust and less powerful.
  Default: `FALSE`.

- FDR:

  Numeric value with the adjusted p-value cut-off reported in the
  output. Default: `0.05`.

- adjustMethod:

  String with the multiple testing correction across the sets. Default:
  `"BH"`.

- regionSets:

  Character vector with the names of the sets to test. Default: `NULL`,
  all of them.

- carryCounts:

  Logical value to indicate whether the counts must travel inside the
  result, so that
  [`plotSetSignal`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetSignal.md)
  can draw the signal without being handed the counts object again.
  Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.setResults` object, or a `RegionSetDE.setResultsList`
when `contrast` is a named list.

## Details

The effect size, not the p-value, is the primary output here. A set of
30,000 promoters tested as if its regions were independent returns a
p-value below anything a computer will print for a mean shift of 0.05
log2, which says nothing about whether the shift matters. The regions of
a set are not independent either, since neighbouring elements inside the
same domain move together, so the variance of the mean log2 fold change
is inflated by the factor `1 + (n - 1) * rho`, with `rho` estimated from
the residuals of the fit through
[`limma::interGeneCorrelation`](https://rdrr.io/pkg/limma/man/camera.html).
The confidence interval in the output carries that inflation; read it
before reading the p-value. The two tests answer different questions and
the pattern between them is informative. `camera` is competitive: it
asks whether the regions of the set moved more than the regions they are
compared against, and it is invariant to a scaling error affecting every
region equally. `fry` is self-contained: it asks whether they moved away
from zero at all, which a global shift in the mark, or a residual
normalisation error, will satisfy for every set at once. When camera
separates the sets and fry does not, the sets redistributed the signal
between them; when fry is significant everywhere and camera nowhere,
everything moved together and the normalisation deserves a second look
before the biology does. The comparison universe comes from the fit,
which built it once, and travels on into the result, so
[`plotUniverseMatching`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotUniverseMatching.md)
can check the matching afterwards without anything being kept on the
side. Passing a `RegionSetDE.universe` object, or one of the two
keywords, overrides it for this test alone. A fit with no replicates
loses the self-contained test. `fry` builds a linear model inside each
set and needs a residual to measure it against, which a design with one
sample per level does not have, so it is dropped with a message and only
the competitive test runs. The correlation between regions goes the same
way: it is estimated from the residuals of the fit, and without them it
falls back to 0.01, the value `limma` uses when nothing better is
available. That number sets how much the confidence interval is widened,
so on such a fit the interval is as assumed as the dispersion is, and
`interRegionCor` is worth setting by hand from a replicated experiment
on the same assay when one exists. The competitive test runs through
[`limma::cameraPR`](https://rdrr.io/pkg/limma/man/camera.html) on the
per-region statistics, which is what makes it work identically for the
four engines. The self-contained test needs the values themselves and is
computed on the log-CPM matrix of the fit; for `edgeR` and `DESeq2` that
matrix is a transformation of the counts rather than the quantity the
model was fitted on, so the two are close but not identical, and the
competitive test is the one to lead with.

## See also

[`testSetContrast`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md),
[`makeSetUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md),
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
[`plotSetEffect`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetEffect.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

setResults <- testRegionSets(fit, contrast = c("condition", "SHR", "BN"),
                             verbose = FALSE)
setResults
#> An object of class 'RegionSetDE.setResults'
#>   test            : set 
#>   contrast        : condition: SHR vs BN 
#>   engine          : edgeR 
#>   methods         : camera, fry 
#>   universe        : otherSets (matched on width and abundance) 
#> 
#>      region.set n.regions mean.log2FC delta.log2FC CI.lower CI.upper camera.FDR
#>     promoterCpG       269     -0.8340      -0.9310   -1.870  0.00836      0.938
#>        geneBody       909      0.2910       0.4070   -0.832  1.64000      0.938
#>      intergenic       440      0.2650       0.2060   -1.140  1.55000      0.938
#>  promoterNonCpG       277     -0.0209      -0.0787   -1.510  1.35000      0.938
#>  fry.FDR
#>     0.54
#>     0.54
#>     0.54
#>     0.54

# The intergenic set is the control: it should not come out as responding
resultsTable(setResults)
#>       region.set n.regions n.comparison mean.log2FC median.log2FC
#> 1    promoterCpG       269          314 -0.83367077    -0.9443600
#> 2       geneBody       909          986  0.29148227     0.1692341
#> 3     intergenic       440         1354  0.26499294     0.1740960
#> 4 promoterNonCpG       277         1378 -0.02089479    -0.1571407
#>   mean.log2FC.comparison delta.log2FC   CI.lower    CI.upper inter.region.cor
#> 1             0.09704300  -0.93071378 -1.8697826 0.008355036        0.8583680
#> 2            -0.11505923   0.40654151 -0.8318836 1.644966633        0.2721337
#> 3             0.05862768   0.20636526 -1.1374078 1.550138278        0.2733108
#> 4             0.05781212  -0.07870691 -1.5089372 1.351523422        0.4433140
#>   median.width camera.direction  camera.p fry.direction     fry.p camera.FDR
#> 1         1000             Down 0.3213638          Down 0.2686203  0.9384488
#> 2         1000               Up 0.5187249          Down 0.5401741  0.9384488
#> 3         1000               Up 0.7773973          Down 0.4717232  0.9384488
#> 4         1000             Down 0.9384488          Down 0.4613130  0.9384488
#>     fry.FDR
#> 1 0.5401741
#> 2 0.5401741
#> 3 0.5401741
#> 4 0.5401741

```
