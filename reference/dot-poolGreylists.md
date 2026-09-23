# .poolGreylists

Pools the regions flagged in every input, keeping the positions flagged
by at least `minInputs` of them, and records which inputs flagged each
region.

## Usage

``` r
.poolGreylists(flaggedList, minInputs, genomeInfo)
```

## Arguments

- flaggedList:

  Named list of `GRanges`, one per input, each without overlaps within
  itself.

- minInputs:

  Integer with the number of inputs that must flag a position.

- genomeInfo:

  `Seqinfo` of the windows.

## Value

A `GRanges` with the `n.inputs` and `inputs` metadata columns.

## Author

Sebastian Gregoricchio
