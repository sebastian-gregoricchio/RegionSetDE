# .asRegionSet

Lets the counting functions take a single `GRanges` where they expect
region sets. The ranges go through
[`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md)
and come back as a `RegionSetDE` object with one set, so that everything
downstream finds the object it expects. Anything else is returned as it
is.

## Usage

``` r
.asRegionSet(regionSet, setName = "regions", verbose = TRUE)
```

## Arguments

- regionSet:

  Object passed to the counting functions.

- setName:

  String with the name given to the set. Default: `"regions"`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE` object when `regionSet` is a `GRanges`, `regionSet`
itself otherwise.

## Author

Sebastian Gregoricchio
