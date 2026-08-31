# .prepareResultTable

Pulls the table out of a results object and restricts it to a subset of
region sets, relabelling `diff.status` when the thresholds differ from
the stored ones.

## Usage

``` r
.prepareResultTable(results, set = NULL, FDR = NULL, log2FC = NULL)
```

## Arguments

- results:

  `RegionSetDE.results` object.

- set:

  Character vector with the region sets to keep, or `NULL`.

- FDR:

  Numeric value with the adjusted p-value cut-off, or `NULL`.

- log2FC:

  Numeric value with the log2 fold change cut-off, or `NULL`.

## Value

A data.frame ready to be plotted.

## Author

Sebastian Gregoricchio
