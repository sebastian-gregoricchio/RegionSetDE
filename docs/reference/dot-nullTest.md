# .nullTest

Fits and tests the null rows the same way the fit tested the regions.

## Usage

``` r
.nullTest(
  dgeList,
  design,
  contrastVector,
  dispersion,
  useQuasiLikelihood = FALSE
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

  Numeric value with the dispersion of the fit.

- useQuasiLikelihood:

  Logical value indicating whether the fit used the quasi-likelihood F
  test.

## Value

A data.frame with the `logFC` and `PValue` columns.

## Author

Sebastian Gregoricchio
