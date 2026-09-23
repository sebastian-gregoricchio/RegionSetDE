# .tileSummaryRule

Decides how the tiles of a region are combined, from the argument when
it is given and from the way the object was counted otherwise.

## Usage

``` r
.tileSummaryRule(counts, tileSummary = NULL)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- tileSummary:

  String with the rule requested, or `NULL`.

## Value

A string, one among `"sum"`, `"mean"`, `"max"` and `"min"`.

## Author

Sebastian Gregoricchio
