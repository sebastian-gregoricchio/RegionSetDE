# .groupAverageSignal

Computes the average signal of every region over the samples of each
level of a column of the `colData`, with the same function the engine
uses for `average.signal` over all the samples, so that the columns can
be read side by side:
[`edgeR::aveLogCPM`](https://rdrr.io/pkg/edgeR/man/aveLogCPM.html) for
edgeR, the mean of the log2 values of the linear model for voom, limma
and dream, and `log2` of the mean normalised count plus one for DESeq2.

## Usage

``` r
.groupAverageSignal(fit, column)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

- column:

  String with the name of the column, or `NULL`.

## Value

A data.frame with one `average.signal.<level>` column per level and one
row per row of the fit, or a data.frame with no column when `column` is
`NULL`.

## Author

Sebastian Gregoricchio
