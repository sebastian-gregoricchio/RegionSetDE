# .testDESeq2

Runs the Wald test on a `DESeq2` fit.

## Usage

``` r
.testDESeq2(fit, contrastVector, lfcThreshold = 0)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

- contrastVector:

  Numeric vector with the contrast.

- lfcThreshold:

  Numeric value with the log2 fold change of the null hypothesis.

## Value

A data.frame with the `log2FC`, `average.signal`, `stat` and `p.value`
columns.

## Author

Sebastian Gregoricchio
