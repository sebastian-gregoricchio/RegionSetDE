# .expressionMatrix

Returns the matrix of log2 values the fit was built on, or the closest
transformation of the counts when the engine works on the count scale.

## Usage

``` r
.expressionMatrix(fit, priorCount = 2)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

- priorCount:

  Numeric value with the prior count added before taking the logarithm.
  Default: `2`.

## Value

A numeric matrix with one row per region and one column per sample.

## Author

Sebastian Gregoricchio
