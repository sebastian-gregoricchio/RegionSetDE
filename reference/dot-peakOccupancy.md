# .peakOccupancy

Tells, for every region, which groups have a consensus peak on it and
how many samples have a peak on it.

## Usage

``` r
.peakOccupancy(regionRanges, groupConsensus, peakList)
```

## Arguments

- regionRanges:

  `GRanges` with the regions.

- groupConsensus:

  Named list of `GRanges`, the consensus of every group.

- peakList:

  Named `GRangesList` with the peaks of every sample.

## Value

A `DataFrame` with one logical `peak.<group>` column per group,
`peak.groups` and `peak.samples`.

## Author

Sebastian Gregoricchio
