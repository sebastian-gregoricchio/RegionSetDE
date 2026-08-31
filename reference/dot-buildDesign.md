# .buildDesign

Turns the `design` argument of
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
into a model matrix, keeping the formula aside when the engine needs it.

## Usage

``` r
.buildDesign(design, random = NULL, colTable, engine)
```

## Arguments

- design:

  Formula, string holding a formula, or design matrix.

- random:

  Formula, or string holding a formula, with the random terms, or
  `NULL`.

- colTable:

  Data.frame with the sample metadata.

- engine:

  String with the engine being used.

## Value

A list with the `matrix` of the fixed effects and the `formula`, `NULL`
when a matrix was supplied.

## Author

Sebastian Gregoricchio
