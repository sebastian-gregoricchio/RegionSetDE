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
  universeSets = NULL,
  effectMethod = "sample",
  interRegionCor = NULL,
  tileHandling = "collapse",
  overlapPolicy = "drop",
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

- universeSets:

  Character vector with the names of the sets the comparison rows are
  drawn from, when one has to be built here. Default: `NULL`, every set
  other than the one being tested.

- effectMethod:

  String with what the confidence interval on the effect is computed
  from, either `"sample"`, which treats the biological samples as the
  replication, or `"region"`, which treats the regions as it. Default:
  `"sample"`.

- interRegionCor:

  Numeric value with the correlation between regions, used to inflate
  the variance of the region-heterogeneity interval of both the set and
  the rows it is compared against. Default: `NULL`, estimated separately
  for each of the two from the residuals of the fit, or held at 0.01
  when the design leaves no residual to estimate it from.

- tileHandling:

  String with what to do when the fit was built on tiles, either
  `"collapse"`, which averages the tiles of a region back into one row
  before the set is assembled, or `"keep"`, which lets every tile count
  on its own. Default: `"collapse"`.

- overlapPolicy:

  String with what to do about comparison rows overlapping the set in
  the genome, one among `"allow"`, `"drop"` and `"stop"`. Default:
  `"drop"`.

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
when `contrast` is a named list. The table carries `CI.lower` and
`CI.upper` for the interval selected by `effectMethod`, `CI.type` naming
which one that is, and `heterogeneity.CI.lower` and
`heterogeneity.CI.upper` for the region-level one, always.

## Details

The effect size, not the p-value, is the primary output here. A set of
30,000 promoters tested as if its regions were independent returns a
p-value below anything a computer will print for a mean shift of 0.05
log2, which says nothing about whether the shift matters.

Two different intervals can be put around that effect and they answer
different questions, so both are reported and `effectMethod` decides
which one is called the confidence interval. The `"sample"` interval is
the default and is the one to quote as a biological result. One number
is computed per library, the mean signal over the set minus the mean
signal over its comparison, and those numbers are then run through the
design of the experiment. The replication is the biological samples,
which is where it comes from in the experiment, and adding regions to a
set makes that interval more stable without ever making it narrower than
four libraries can support.

The `"region"` interval is the mean of the per-region log2 fold changes
with its variance inflated by `1 + (n - 1) * rho`, with `rho` estimated
from the residuals of the fit through
[`limma::interGeneCorrelation`](https://rdrr.io/pkg/limma/man/camera.html).
It describes how much the effect varies from locus to locus within the
set, conditional on these libraries, and it is reported under
`heterogeneity.CI.lower` and `heterogeneity.CI.upper` whatever
`effectMethod` is set to. It is a useful quantity and it is not a
confidence interval on a condition effect: the sampling units behind it
are genomic loci, and no number of loci substitutes for the biological
replication that was or was not done. Reading it as the second thing
rather than the first is the safe habit. Note also that its width does
not fall away as the set grows, since `sd^2 / n * (1 + (n - 1) * rho)`
tends to `sd^2 * rho`; it flattens rather than collapsing.

The two tests answer different questions and neither of them, alone or
in combination, establishes that a set did not change. `camera` is
competitive: it asks whether the regions of the set moved more than the
regions they are compared against, and it is invariant to a scaling
error affecting every region equally. `fry` is self-contained: it asks
whether they moved away from zero at all. Read the four outcomes as
evidence and not as mechanism:

- camera significant, fry significant: evidence both that the set moved
  and that it moved more than its comparison.

- camera significant, fry not: evidence of a difference relative to the
  comparison, with the absolute claim left open. Failing to reject the
  self-contained null is not evidence that the absolute change is zero,
  and the two tests do not have the same power.

- camera not significant, fry significant: evidence that the set moved,
  none that it moved differently from its comparison.

- neither significant: neither test found evidence, which is not the
  same as evidence of no effect.

The word for a set gaining what another set lost is redistribution, and
it is a claim about two sets rather than about one set and its universe,
so it belongs to
[`testSetContrast`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md).
To argue that a mark did not change globally, an equivalence test
against a bounded near-zero interval is what the claim needs; a
non-significant `fry` is not that. Bear in mind too that any centring
normalisation, TMM and background included, removes a genuinely global
shift from the data before `fry` ever sees it, so the self-contained
test is not the place to look for one.

The comparison universe comes from the fit, which built it once, and
travels on into the result, so
[`plotUniverseMatching`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotUniverseMatching.md)
can check the matching afterwards without anything being kept on the
side. Passing a `RegionSetDE.universe` object, or one of the two
keywords, overrides it for this test alone. Whatever it is, it is made
of the other sets loaded into the object, and the competitive p-value is
a statement about the set relative to those and not relative to the
genome. Load two sets and the test compares them to each other; load
four that behave alike and every one of them can come out unremarkable
against the other three.
[`makeSetUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md)
takes `universeSets` for choosing that comparison pool explicitly, which
is worth doing when the sets were not all picked for the same reason.

Both intervals carry the uncertainty of both sides. The regions of the
comparison are no less correlated than the regions of the set, so
treating their mean as if it were known would leave the interval
narrower than the data supports.

A set that overlaps its own comparison in the genome shares reads with
it and drags the difference towards zero. Overlap is measured on the
coordinates rather than on the identifiers, so two sets holding
chr1:1000-2000 and chr1:1500-2500 are seen as overlapping even though no
region identifier is shared, and `overlapPolicy` decides what happens
next. The number of comparison rows removed, or left in place, is
reported in `n.comparison.overlapping`.

On a tiled fit the row is a tile, and a set assembled from tiles weights
each region by how many tiles it was cut into: a 40 kb domain would
count forty times a 2 kb one. That changes the question from the average
response of the regions in the set to the average response of the base
pairs in it. `tileHandling = "collapse"` averages the tiles of a region
back together first, which keeps the region as the unit and matches what
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
does at its own level. `"keep"` is the base-pair version, and is a
deliberate choice rather than a default.

A fit with no replicates loses the self-contained test. `fry` builds a
linear model inside each set and needs a residual to measure it against,
which a design with one sample per level does not have, so it is dropped
with a message and only the competitive test runs. The correlation
between regions goes the same way: it is estimated from the residuals of
the fit, and without them it falls back to 0.01, the value `limma` uses
when nothing better is available. That number sets how much the
confidence interval is widened, so on such a fit the interval is as
assumed as the dispersion is, and `interRegionCor` is worth setting by
hand from a replicated experiment on the same assay when one exists. The
competitive test runs through
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

# The universe comes from the fit and travels into the result
setRes <- testRegionSets(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
resultsTable(setRes)
#>       region.set n.regions n.comparison n.comparison.overlapping mean.log2FC
#> 1    promoterCpG       269          314                        0 -0.83367077
#> 2       geneBody       909          986                        0  0.29148227
#> 3     intergenic       440         1354                        0  0.26499294
#> 4 promoterNonCpG       277         1378                        0 -0.02089479
#>   median.log2FC mean.log2FC.comparison delta.log2FC   CI.lower  CI.upper
#> 1    -0.9443600             0.09704300  -0.93071378 -1.9314597 0.5607824
#> 2     0.1692341            -0.11505923   0.40654151 -0.3328615 0.9187083
#> 3     0.1740960             0.05862768   0.20636526 -0.1710428 0.3668162
#> 4    -0.1571407             0.05781212  -0.07870691 -0.2633465 0.2509689
#>   CI.type heterogeneity.CI.lower heterogeneity.CI.upper inter.region.cor
#> 1  sample              -2.547867              0.6864395        0.8583680
#> 2  sample              -1.574041              2.3871240        0.2721337
#> 3  sample              -1.761069              2.1737990        0.2733108
#> 4  sample              -2.067061              1.9096471        0.4433140
#>   inter.region.cor.universe median.width sample.delta.log2FC sample.delta.SE
#> 1                 0.3775802         1000        -0.685338676      0.28961692
#> 2                 0.4530172         1000         0.292923395      0.14544165
#> 3                 0.3815880         1000         0.097886686      0.06250318
#> 4                 0.3504707         1000        -0.006188815      0.05976724
#>   sample.delta.df sample.delta.p camera.direction  camera.p fry.direction
#> 1               2      0.1416115             Down 0.3213638          Down
#> 2               2      0.1816079               Up 0.5187249          Down
#> 3               2      0.2578184               Up 0.7773973          Down
#> 4               2      0.9269756             Down 0.9384488          Down
#>       fry.p camera.FDR   fry.FDR sample.delta.FDR
#> 1 0.2686203  0.9384488 0.5401741        0.3437579
#> 2 0.5401741  0.9384488 0.5401741        0.3437579
#> 3 0.4717232  0.9384488 0.5401741        0.3437579
#> 4 0.4613130  0.9384488 0.5401741        0.9269756

plotUniverseMatching(setRes)

plotSetEffect(setRes)


# Overriding it for one test
setRes <- testRegionSets(fit, contrast = c("condition", "SHR", "BN"), universe = "all", verbose = FALSE)
```
