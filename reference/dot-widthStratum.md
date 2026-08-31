# .widthStratum

Splits a vector of widths into strata of comparable size, used to keep
the relative filters from favouring the wide regions.

## Usage

``` r
.widthStratum(rowWidths, strataNumber = 5)
```

## Arguments

- rowWidths:

  Numeric vector with the widths.

- strataNumber:

  Numeric value with the number of strata.

## Value

A character vector with the stratum of every row.

## Author

Sebastian Gregoricchio
