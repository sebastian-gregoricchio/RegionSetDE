# .setStatistics

Runs the per-region contrast and returns the statistics the set level
tests are built on, one row per row of the fit.

## Usage

``` r
.setStatistics(fit, contrastObject)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

- contrastObject:

  List returned by `.resolveContrast`.

## Value

A data.frame with the region annotation and the per-region `log2FC` and
`stat` columns.

## Author

Sebastian Gregoricchio
