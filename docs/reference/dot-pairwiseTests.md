# .pairwiseTests

Runs a two-sample test between each pair of groups on the values of one
region, as they are drawn, pairing the samples through a column of the
metadata when asked to.

## Usage

``` r
.pairwiseTests(
  plotTable,
  colTable,
  test = "t.test",
  pairBy = NULL,
  comparisons = NULL,
  groupLevels,
  pAdjustMethod = "none"
)
```

## Arguments

- plotTable:

  Data.frame with one row per sample and the `sample`, `value` and
  `group` columns.

- colTable:

  Data.frame with the sample metadata and a `sample` column.

- test:

  String with the test, either `"t.test"` or `"wilcox.test"`. Default:
  `"t.test"`.

- pairBy:

  String with the column pairing the samples, or `NULL` for unpaired
  tests. Default: `NULL`.

- comparisons:

  List of character vectors of length two, or `NULL` for every pair of
  groups. Default: `NULL`.

- groupLevels:

  Character vector with the groups, in the order of the axis.

- pAdjustMethod:

  String with the correction applied across the comparisons. Default:
  `"none"`.

## Value

A data.frame with one row per comparison and the `group1`, `group2`,
`n.first`, `n.second`, `p.value` and `p.adjusted` columns. With `pairBy`
the two sizes are the number of pairs.

## Author

Sebastian Gregoricchio
