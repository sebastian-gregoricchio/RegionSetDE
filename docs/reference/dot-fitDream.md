# .fitDream

Fits the mixed model of `variancePartition`, with precision weights
estimated under the same formula.

## Usage

``` r
.fitDream(
  countMatrix,
  designMatrix,
  designFormula,
  librarySizes,
  offsetMatrix,
  colTable,
  BPPARAM = NULL
)
```

## Arguments

- countMatrix:

  Numeric matrix of counts.

- designMatrix:

  Design matrix of the fixed effects, kept for the contrast resolution.

- designFormula:

  Formula including the random terms.

- librarySizes:

  Numeric vector with the library sizes.

- offsetMatrix:

  Matrix of log offsets, or `NULL`.

- colTable:

  Data.frame with the sample metadata.

- BPPARAM:

  `BiocParallelParam` object, or `NULL`.

## Value

A list with the `fit`, the `blocking` and the `dispersion` elements.

## Author

Sebastian Gregoricchio
