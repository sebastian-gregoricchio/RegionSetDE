# .interRegionCor

Estimates the correlation between the regions of a set, from the
residuals of the design.

## Usage

``` r
.interRegionCor(expressionMatrix, design, index, label = "the set")
```

## Arguments

- expressionMatrix:

  Numeric matrix of log2 values.

- design:

  Design matrix.

- index:

  Integer vector with the rows of the set.

- label:

  String naming what is being estimated, used in the warning. Default:
  `"the set"`.

## Value

A numeric value.

## Details

When the design leaves fewer than two residual degrees of freedom there
is nothing to estimate a correlation from, and the value falls back to
the 0.01 that `limma` uses in the same situation. That fallback is loud
rather than silent, because the difference between 0.01 and a measured
0.4 is a fortyfold change in every variance, and a design that gains a
coefficient can cross that line without anything else about the analysis
appearing to change.

## Author

Sebastian Gregoricchio
