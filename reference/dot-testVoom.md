# .testVoom

Runs the moderated t test on a `limma` fit.

## Usage

``` r
.testVoom(fit, contrastVector, lfcThreshold = 0, trend = FALSE)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

- contrastVector:

  Numeric vector with the contrast.

- lfcThreshold:

  Numeric value with the log2 fold change of the null hypothesis.

## Value

A data.frame with the `log2FC`, `average.signal`, `stat`,
`stat.distribution`, `df1`, `df2` and `p.value` columns.

## Author

Sebastian Gregoricchio
