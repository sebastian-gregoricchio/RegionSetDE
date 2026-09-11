# .fitLimma

Fits the `limma` model on the log2 signal, with the residual variance
trended on abundance, for values that are not counts.

## Usage

``` r
.fitLimma(
  countMatrix,
  designMatrix,
  librarySizes,
  offsetMatrix,
  colTable,
  block = NULL,
  robust = TRUE,
  priorCount = 2,
  verbose = TRUE
)
```

## Arguments

- countMatrix:

  Numeric matrix with the values.

- designMatrix:

  Design matrix.

- librarySizes:

  Numeric vector with the library size of each sample.

- offsetMatrix:

  Numeric matrix of log offsets, or `NULL`.

- colTable:

  Data.frame with the sample annotation.

- block:

  String naming a `colData` column holding a blocking variable, or
  `NULL`. Default: `NULL`.

- robust:

  Logical value to indicate whether the prior variance must be estimated
  robustly. Default: `TRUE`.

- priorCount:

  Numeric value added before the logarithm. Default: `2`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A list with the fit, the blocking information and the median residual
standard deviation.

## Details

Where `voom` estimates a weight per observation from a mean-variance
trend fitted on the count scale, this fits the trend on the residual
variance itself and applies no weights. Nothing about the values then
has to be a fragment count, which is what makes it the engine for
coverage, for an externally normalised matrix, or for anything else that
arrives already on a continuous scale.

## Author

Sebastian Gregoricchio
