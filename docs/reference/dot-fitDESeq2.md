# .fitDESeq2

Fits the negative binomial model of `DESeq2` on a user supplied model
matrix.

## Usage

``` r
.fitDESeq2(
  countMatrix,
  designMatrix,
  offsetMatrix,
  librarySizes,
  colTable,
  verbose = TRUE
)
```

## Arguments

- countMatrix:

  Integer matrix of counts.

- designMatrix:

  Design matrix.

- offsetMatrix:

  Matrix of log offsets, or `NULL`.

- librarySizes:

  Numeric vector with the library sizes.

- colTable:

  Data.frame with the sample metadata.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A list with the `fit`, the `blocking` and the `dispersion` elements.

## Author

Sebastian Gregoricchio
