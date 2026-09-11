# .testDream

Runs the test on a `dream` fit. A contrast that is not a single
coefficient of the design needs the mixed model to be fitted again,
since `dream` builds the contrast at fit time.

## Usage

``` r
.testDream(fit, contrastObject, lfcThreshold = 0, verbose = TRUE)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

- contrastObject:

  List returned by `.resolveContrast`.

- lfcThreshold:

  Numeric value with the log2 fold change of the null hypothesis.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A data.frame with the `log2FC`, `average.signal`, `stat` and `p.value`
columns.

## Author

Sebastian Gregoricchio
