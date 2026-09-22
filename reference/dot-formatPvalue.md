# .formatPvalue

Formats p-values for a markdown label, with a fixed number of decimals
from 0.1 upwards and in scientific notation below, e.g. 3.20 x 10 to the
-2 written with a superscript exponent. A value that underflowed to zero
is written as a bound.

## Usage

``` r
.formatPvalue(p, decimals = 2)
```

## Arguments

- p:

  Numeric vector with the p-values.

- decimals:

  Numeric value with the number of decimals. Default: `2`.

## Value

A character vector with one formatted value per p-value, `"NA"` for the
missing ones.

## Author

Sebastian Gregoricchio
