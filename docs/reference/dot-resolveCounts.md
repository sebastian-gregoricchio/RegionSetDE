# .resolveCounts

Returns the counts and, when there are any, the statistics a plotting
function has to work with, whatever kind of object it was handed.

## Usage

``` r
.resolveCounts(object, counts = NULL, contrast = NULL)
```

## Arguments

- object:

  Any object of the package holding counts: `RegionSetDE.counts`,
  `RegionSetDE.fit`, `RegionSetDE.results`, `RegionSetDE.setResults`, or
  one of the two list classes.

- counts:

  `RegionSetDE.counts` object overriding the one carried by `object`.
  Default: `NULL`.

- contrast:

  String with the name of a contrast, or its position. Default: `NULL`.

## Value

A list with the `counts` and the `results` elements, the second one
`NULL` when no statistics are available.

## Author

Sebastian Gregoricchio
