# .resolveUniverse

Turns the `universe` argument of
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
and
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
into a `RegionSetDE.universe` object, building it when a keyword was
given rather than an object.

## Usage

``` r
.resolveUniverse(
  object,
  universe = "matched",
  matchOn = c("width", "abundance"),
  universeRatio = 5,
  regionSets = NULL,
  soft = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  `RegionSetDE.fit` or `RegionSetDE.counts` object.

- universe:

  String with a keyword, a `RegionSetDE.universe` object, or `NULL`.

- matchOn:

  Character vector with the covariates the comparison rows are matched
  on.

- universeRatio:

  Numeric value with the number of comparison rows drawn per region of
  the set.

- regionSets:

  Character vector with the sets being tested, or `NULL`.

- soft:

  Logical value to indicate whether a universe that cannot be built must
  return empty rather than stop. Default: `FALSE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.universe` object, empty when none could be built under
`soft`.

## Author

Sebastian Gregoricchio
