# .modelBrackets

Reads, for every contrast comparing two levels of the grouping column,
the log2 fold change and the FDR the fit gave to one region, which is
what the model brackets of
[`plotRegion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegion.md)
write.

## Usage

``` r
.modelBrackets(
  resultsList,
  regionKey,
  groupBy,
  groupLevels,
  comparisons = NULL
)
```

## Arguments

- resultsList:

  Named list of `RegionSetDE.results` objects.

- regionKey:

  String with the region, written as `"set|id"`.

- groupBy:

  String with the column grouping the samples on the axis.

- groupLevels:

  Character vector with the groups, in the order of the axis.

- comparisons:

  List of character vectors of length two, or `NULL`. Default: `NULL`.

## Value

A data.frame with one row per bracket and the `contrast`, `group1`,
`group2`, `log2FC` and `FDR` columns.

## Author

Sebastian Gregoricchio
