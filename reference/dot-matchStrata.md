# .matchStrata

Draws comparison rows from the strata occupied by a region set, so that
the two share a distribution of width and abundance.

## Usage

``` r
.matchStrata(
  setRows,
  eligibleRows,
  match = c("width", "abundance"),
  ratio = 5,
  strata = 5
)
```

## Arguments

- setRows:

  Data.frame with the rows of the set.

- eligibleRows:

  Data.frame with the rows the comparison can be drawn from.

- match:

  Character vector with the covariates to match on.

- ratio:

  Numeric value with the number of comparison rows drawn per region of
  the set.

- strata:

  Numeric value with the number of strata per covariate.

## Value

An integer vector with the positions of the comparison rows.

## Author

Sebastian Gregoricchio
