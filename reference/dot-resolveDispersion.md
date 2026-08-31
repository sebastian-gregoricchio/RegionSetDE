# .resolveDispersion

Works out which dispersion
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
should hold fixed, estimating one from rows assumed not to respond when
the design leaves no residual to take it from.

## Usage

``` r
.resolveDispersion(
  counts,
  dispersion = NULL,
  residualDegrees,
  engine = "edgeR",
  nullSource = "background",
  nullRegionSets = NULL,
  verbose = TRUE
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- dispersion:

  Numeric value, list from
  [`estimateNullDispersion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md),
  keyword, or `NULL`.

- residualDegrees:

  Numeric value with the residual degrees of freedom of the design.

- engine:

  String with the engine being used.

- nullSource:

  String with where the null rows come from.

- nullRegionSets:

  Character vector with the sets used as null rows, or `NULL`.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A list with the `dispersion`, `NULL` when it must come from the residual
variation, the `source` it came from, and the rows held out of the
estimate.

## Author

Sebastian Gregoricchio
