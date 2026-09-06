# .overlappingRows

Finds the comparison rows that overlap a region set in the genome.

## Usage

``` r
.overlappingRows(regionStats, setIndex, comparisonIndex)
```

## Arguments

- regionStats:

  Data.frame returned by `.setStatistics`.

- setIndex:

  Integer vector with the rows of the set.

- comparisonIndex:

  Integer vector with the rows the set is compared against.

## Value

An integer vector with the positions of the overlapping comparison rows.

## Details

Two region sets can describe the same chromatin without sharing a single
identifier, so the check is on the coordinates. A comparison row that
overlaps the set carries some of the same reads as the set does, which
pulls the difference between the two towards zero however the
identifiers were assigned.

## Author

Sebastian Gregoricchio
