# .fitVoom

Fits the `limma` model on log-CPM values with the precision weights of
`voom`.

## Usage

``` r
.fitVoom(
  countMatrix,
  designMatrix,
  librarySizes,
  offsetMatrix,
  colTable,
  block = NULL,
  robust = TRUE,
  verbose = TRUE
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

- block:

  String with the name of the blocking column, or `NULL`.

- robust:

  Logical value passed to the empirical Bayes step.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A list with the `fit`, the `blocking` and the `dispersion` elements.

## Author

Sebastian Gregoricchio
