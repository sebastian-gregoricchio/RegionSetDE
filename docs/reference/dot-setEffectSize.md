# .setEffectSize

Computes the mean log2 fold change of a set, the difference with its
background, and a region-heterogeneity interval inflated for the
correlation between the regions.

## Usage

``` r
.setEffectSize(
  logFC,
  setIndex,
  backgroundIndex,
  correlation,
  backgroundCorrelation = NULL,
  level = 0.95
)
```

## Arguments

- logFC:

  Numeric vector with the per-region log2 fold changes.

- setIndex:

  Integer vector with the rows of the set.

- backgroundIndex:

  Integer vector with the rows of the background.

- correlation:

  Numeric value with the correlation between the regions of the set.

- backgroundCorrelation:

  Numeric value with the correlation between the regions of the
  background. Default: `NULL`, the same as the set, since a comparison
  drawn from the same object is correlated in the same way.

- level:

  Numeric value with the confidence level. Default: `0.95`.

## Value

A list with the means, the difference and the bounds of the interval.

## Details

The sampling units here are the genomic loci, so the interval describes
how much the effect varies from one region of the set to another. It is
not a confidence interval on a condition effect, whatever the number of
regions: the biological replication of the experiment lives in the
samples, and `.sampleSetEffect` is where that interval comes from. The
two are reported side by side and `effectMethod` decides which one
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
labels as the confidence interval.

Growing the set does not shrink this interval away. With `rho` held
constant, `sd^2 / n * (1 + (n - 1) * rho)` tends to `sd^2 * rho` as the
set grows, so the width flattens out at a floor set by the correlation
rather than falling towards zero. What the inflation does not fix is
that `sd` is the spread of estimated fold changes, which carries the
estimation error of each region along with the genuine variation between
them.

## Author

Sebastian Gregoricchio
