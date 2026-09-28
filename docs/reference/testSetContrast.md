# testSetContrast

Asks whether a contrast affects one region set differently from another.
This is the comparison behind questions of the kind "does the treatment
reduce the mark more at Polycomb promoters than at active enhancers",
and it is the one claim that a global normalisation error cannot
manufacture, since a scaling factor that is wrong for one set is wrong
for the other in the same way.

## Usage

``` r
testSetContrast(
  fit,
  contrast,
  set1 = NULL,
  set2 = NULL,
  effectMethod = "sample",
  interRegionCor = NULL,
  useRanks = FALSE,
  sharedRegions = "drop",
  overlapWithinSet = "warn",
  FDR = 0.05,
  adjustMethod = "BH",
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

- set1:

  Character vector with the name, or names, of the first region set.
  Default: `NULL`, every pair of sets is tested.

- set2:

  Character vector with the name, or names, of the second region set.
  Default: `NULL`.

- effectMethod:

  String with what the confidence interval on the difference is computed
  from, either `"sample"` or `"region"`. Default: `"sample"`.

- interRegionCor:

  Numeric value with the correlation between the regions of a set.
  Default: `NULL`, estimated from the residuals.

- useRanks:

  Logical value to indicate whether the test must work on the ranks
  rather than on the statistics. Default: `FALSE`.

- sharedRegions:

  String with what to do with the regions the two sets share in the
  genome, either `"drop"` or `"stop"`. Default: `"drop"`.

- overlapWithinSet:

  String with what to do when regions on the same side of the pair
  overlap each other in the genome, one among `"allow"`, `"warn"` and
  `"stop"`. The counts are reported in `n.overlapping.within.1` and
  `n.overlapping.within.2`. Default: `"warn"`.

- FDR:

  Numeric value with the adjusted p-value cut-off reported in the
  output. Default: `0.05`.

- adjustMethod:

  String with the multiple testing correction across the pairs. Default:
  `"BH"`.

- carryCounts:

  Logical value to indicate whether the counts must travel inside the
  result. Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.setResults` object with one row per pair of sets, or a
`RegionSetDE.setResultsList` when `contrast` is a named list.

## Details

The test restricts the universe to the two sets and runs the competitive
test of
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
on the first of them, which is exactly a comparison of the first set
against the second. The effect size is the difference between the two
mean log2 fold changes, with the interval selected by `effectMethod`
beside it and the region-heterogeneity one always reported next to it.
As in
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md),
`mean.log2FC.1` and `mean.log2FC.2` are means of per-region fold
changes, every region counting once, and not the fold changes of the
signal pooled over each set.

This is the function for the redistribution question. A set gaining what
another set lost is a claim about two sets, and it is the one claim a
global normalisation error cannot manufacture, since a scaling factor
that is wrong for one set is wrong for the other in the same way. Asking
it through the pattern of a competitive and a self-contained test on a
single set does not work, because failing to reject a self-contained
null is not evidence that the absolute change was zero.

A region shared by the two sets carries the same reads into both sides
of the comparison and pulls the difference towards zero. Sharing is
measured on the genome and not on the identifiers: chr1:1000-2000 in one
set and chr1:1500-2500 in the other are half the same chromatin even
though neither region identifier appears twice. Overlapping regions are
removed from both sides by default and the number removed is reported in
`n.shared.dropped`; `sharedRegions = "stop"` refuses to run instead,
which is the safer setting when the overlap is unexpected.

Overlaps within one side of the pair are a separate matter, measured
after the shared regions are removed and handled by `overlapWithinSet`
as
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
handles them.

## See also

[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md),
[`plotSetEffect`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetEffect.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

setContrast <- testSetContrast(fit, contrast = c("condition", "SHR", "BN"),
                               set1 = "promoterCpG", set2 = "promoterNonCpG", verbose = FALSE)
resultsTable(setContrast)
#>         set.1          set.2 n.regions.1 n.regions.2 n.shared.dropped
#> 1 promoterCpG promoterNonCpG         269         277                0
#>   n.overlapping.within.1 n.overlapping.within.2 mean.log2FC.1 mean.log2FC.2
#> 1                      0                      0    -0.8336708   -0.02089479
#>   delta.log2FC  CI.lower  CI.upper CI.type heterogeneity.CI.lower
#> 1    -0.812776 -1.706428 0.4466904  sample              -2.518293
#>   heterogeneity.CI.upper inter.region.cor.1 inter.region.cor.2 camera.direction
#> 1              0.8927412           0.858368           0.443314             Down
#>   camera.p sample.delta.log2FC sample.delta.SE sample.delta.df sample.delta.p
#> 1 0.393662           -0.629869       0.2502083               2      0.1281565
#>   camera.FDR sample.delta.FDR
#> 1   0.393662        0.1281565

# Every pair at once
allPairs <- testSetContrast(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
```
