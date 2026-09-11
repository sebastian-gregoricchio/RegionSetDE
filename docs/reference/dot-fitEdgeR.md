# .fitEdgeR

Fits the quasi-likelihood negative binomial model of `edgeR`.

## Usage

``` r
.fitEdgeR(
  countMatrix,
  designMatrix,
  librarySizes,
  offsetMatrix,
  colTable,
  robust = TRUE,
  dispersion = NULL
)
```

## Arguments

- countMatrix:

  Numeric matrix of counts.

- designMatrix:

  Design matrix.

- librarySizes:

  Numeric vector with the library sizes.

- offsetMatrix:

  Matrix of log offsets, or `NULL`.

- colTable:

  Data.frame with the sample metadata.

- robust:

  Logical value passed to the dispersion estimation.

- dispersion:

  Numeric value with a dispersion held fixed, or `NULL` to estimate one.

## Value

A list with the `fit`, the `blocking` and the `dispersion` elements.

## Author

Sebastian Gregoricchio
