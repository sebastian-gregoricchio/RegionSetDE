# .setEffectSize

Computes the mean log2 fold change of a set, the difference with its
background, and a confidence interval inflated for the correlation
between the regions.

## Usage

``` r
.setEffectSize(
  logFC,
  setIndex,
  backgroundIndex,
  correlation,
  backgroundCorrelation = 0,
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
  background. Default: `0`.

- level:

  Numeric value with the confidence level. Default: `0.95`.

## Value

A list with the means, the difference and the bounds of the interval.

## Author

Sebastian Gregoricchio
