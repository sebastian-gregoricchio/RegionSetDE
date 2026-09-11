# .correlationSummary

Summarises a correlation matrix as the mean within a group against the
mean between groups, for the label of a panel.

## Usage

``` r
.correlationSummary(correlationMatrix, groupVector = NULL, digits = 3)
```

## Arguments

- correlationMatrix:

  Numeric matrix of correlations.

- groupVector:

  Character vector with the group of every sample, or `NULL`.

- digits:

  Numeric value with the number of decimals written.

## Value

A string, empty when no grouping was given.

## Author

Sebastian Gregoricchio
