# .loadExclusionRegions

Loads one or several exclusion lists and pools them into a single set of
ranges, in the chromosome style of the peaks.

## Usage

``` r
.loadExclusionRegions(
  excludeRegions,
  seqlevelsStyle = "UCSC",
  listLabel = "blacklist"
)
```

## Arguments

- excludeRegions:

  A `GRanges`, a path, a data.frame, or a list of them.

- seqlevelsStyle:

  String with the chromosome naming style, or `NULL`.

- listLabel:

  String with the name of the parameter the list came from, used in the
  messages. Default: `"blacklist"`.

## Value

A `GRanges` without overlaps.

## Author

Sebastian Gregoricchio
