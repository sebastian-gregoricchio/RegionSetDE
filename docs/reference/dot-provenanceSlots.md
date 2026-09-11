# .provenanceSlots

Collects the provenance slots of a `RegionSetDE` object. Regions
arriving as a plain `GRangesList` carry no history, so the empty
defaults are returned instead.

## Usage

``` r
.provenanceSlots(regionSet)
```

## Arguments

- regionSet:

  Object passed to the counting functions.

## Value

A named list with the slots shared by the `RegionSetDE.provenance`
classes.

## Author

Sebastian Gregoricchio
