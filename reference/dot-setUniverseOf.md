# .setUniverseOf

Picks the comparison universe a set level test must use: the one carried
by the fit, or one built here when the caller asked for something else.

## Usage

``` r
.setUniverseOf(
  fit,
  universe = NULL,
  matchOn = c("width", "abundance"),
  universeRatio = 5,
  regionSets = NULL,
  verbose = TRUE
)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

- universe:

  String with a keyword, a `RegionSetDE.universe` object, or `NULL` to
  take the one in the fit.

- matchOn:

  Character vector with the covariates the comparison rows are matched
  on.

- universeRatio:

  Numeric value with the number of comparison rows drawn per region of
  the set.

- regionSets:

  Character vector with the sets being tested, or `NULL`.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A `RegionSetDE.universe` object.

## Author

Sebastian Gregoricchio
