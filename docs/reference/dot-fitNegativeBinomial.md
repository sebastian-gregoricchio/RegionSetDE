# .fitNegativeBinomial

Fits a negative binomial to window counts by maximum likelihood. The
mean is estimated by the sample mean, which is its maximum likelihood
estimate whatever the size, and the size by maximising the likelihood
over its logarithm, computed once per distinct count value.

## Usage

``` r
.fitNegativeBinomial(windowCounts)
```

## Arguments

- windowCounts:

  Integer vector with the count of every window.

## Value

A list with `mu` and `size`, the latter infinite when the counts are not
overdispersed.

## Author

Sebastian Gregoricchio
