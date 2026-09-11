# .logSignalMatrix

Puts a matrix of values on a log2 scale comparable across samples,
dividing out the offsets when there are any and the library sizes
otherwise.

## Usage

``` r
.logSignalMatrix(
  countMatrix,
  librarySizes,
  offsetMatrix = NULL,
  priorCount = 2
)
```

## Arguments

- countMatrix:

  Numeric matrix with the values.

- librarySizes:

  Numeric vector with the library size of each sample.

- offsetMatrix:

  Numeric matrix of log offsets, or `NULL`.

- priorCount:

  Numeric value added before the logarithm, damping the variance of the
  low rows. Default: `2`.

## Value

A numeric matrix of log2 values on a per-million scale.

## Author

Sebastian Gregoricchio
