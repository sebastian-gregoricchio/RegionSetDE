# .resolveScaleLimits

Works out the range of a colour scale, filling in either end from the
data when it was left open, and saying how many values will be drawn at
the ends rather than at their own position.

## Usage

``` r
.resolveScaleLimits(values, limits = NULL, drawnValues = NULL, label = "value")
```

## Arguments

- values:

  Numeric vector the default range is read from.

- limits:

  Numeric vector of length two, either element possibly `NA`, or `NULL`.

- drawnValues:

  Numeric vector with every value that will be drawn, used to count the
  ones falling outside. Default: `NULL`, `values`.

- label:

  String naming the quantity, used in the message. Default: `"value"`.

## Value

A numeric vector of length two.

## Author

Sebastian Gregoricchio
