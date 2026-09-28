# .profileRegions

Picks the regions of a profile and the point each of them is aligned on:
the differential regions of a result, the region sets of a counts
object, or regions given by the user.

## Usage

``` r
.profileRegions(
  counts,
  results,
  regions,
  set,
  direction,
  FDR,
  log2FC,
  maxRegions,
  centre
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- results:

  `RegionSetDE.results` object, or `NULL`.

- regions:

  Regions given by the user, or `NULL`.

- set:

  Character vector with the region sets used, or `NULL`.

- direction:

  String, one among `"both"`, `"up"` and `"down"`.

- FDR:

  Numeric value with the adjusted p-value cut-off, or `NULL`.

- log2FC:

  Numeric value with the absolute log2 fold change cut-off, or `NULL`.

- maxRegions:

  Numeric value with the maximum number of regions per group of rows.

- centre:

  String, either `"summit"` or `"midpoint"`.

## Value

A `GRanges` with one element per row, carrying `row.group`, `region.key`
and `centre.position`.

## Author

Sebastian Gregoricchio
