# .bracketLabels

Writes the labels of the brackets as markdown, either as numbers or as
significance symbols, with the log2 fold change above them when it is
given.

## Usage

``` r
.bracketLabels(
  pValues,
  prefix = "p",
  pLabel = "p.value",
  pDecimals = 2,
  log2FC = NULL
)
```

## Arguments

- pValues:

  Numeric vector with the values to write.

- prefix:

  String naming the value, e.g. `"p"` or `"FDR"`. Default: `"p"`.

- pLabel:

  String, either `"p.value"` or `"stars"`. Default: `"p.value"`.

- pDecimals:

  Numeric value with the number of decimals. Default: `2`.

- log2FC:

  Numeric vector with the fold changes written on a first line, or
  `NULL`. Default: `NULL`.

## Value

A character vector with one label per value.

## Author

Sebastian Gregoricchio
