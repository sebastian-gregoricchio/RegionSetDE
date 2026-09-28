# .withinSetOverlap

Counts the rows of a set that overlap another row of the same set in the
genome.

## Usage

``` r
.withinSetOverlap(regionStats, index)
```

## Arguments

- regionStats:

  Data.frame returned by `.setStatistics`.

- index:

  Integer vector with the rows of the set.

## Value

An integer value.

## Details

Adjacent tiles of one region touch without sharing a base, so a tiled
set kept at the tile level is not reported here.

## Author

Sebastian Gregoricchio
