# .greenlistFactors

Turns a matrix of greenlist counts into one scaling factor per sample.

## Usage

``` r
.greenlistFactors(greenlistMatrix, estimator = "medianRatio")
```

## Arguments

- greenlistMatrix:

  Numeric matrix with one row per greenlist region and one column per
  sample.

- estimator:

  String with the estimator, one among `"medianRatio"`, `"TMM"` and
  `"sum"`.

## Value

A numeric vector with one factor per sample, centred on one.

## Author

Sebastian Gregoricchio
