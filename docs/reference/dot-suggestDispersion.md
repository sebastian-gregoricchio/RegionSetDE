# .suggestDispersion

Searches for the dispersion that would leave the expected fraction of
the null rows below the p-value cut-off.

## Usage

``` r
.suggestDispersion(
  dgeList,
  design,
  contrastVector,
  dispersion,
  useQuasiLikelihood = FALSE,
  target = 0.05,
  steps = 12
)
```

## Arguments

- dgeList:

  `DGEList` with the null rows.

- design:

  Design matrix.

- contrastVector:

  Numeric vector with the contrast.

- dispersion:

  Numeric value with the dispersion the fit used, which anchors the
  search.

- useQuasiLikelihood:

  Logical value indicating whether the fit used the quasi-likelihood F
  test.

- target:

  Numeric value with the proportion aimed at.

- steps:

  Numeric value with the number of dispersions tried. Default: `12`.

## Value

A numeric value, `NA` when the target lies outside the range searched.

## Author

Sebastian Gregoricchio
