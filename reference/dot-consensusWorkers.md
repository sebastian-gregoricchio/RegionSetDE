# .consensusWorkers

Shares the cores between the groups built in parallel and the
calibration run inside each group.

## Usage

``` r
.consensusWorkers(BPPARAM, nGroups, calibrate = FALSE)
```

## Arguments

- BPPARAM:

  A `BiocParallelParam` object, a number of cores, or `NULL` for one
  core.

- nGroups:

  Number of groups needing a consensus, the groups of a single sample
  left out.

- calibrate:

  Logical value to indicate whether the threshold is calibrated inside
  each group.

## Value

A list with `outer`, the back end or the number of cores running the
groups, `inner`, the back end or the number of cores each group gets for
its permutations, and `nOuter` and `nInner`, the number of workers of
each.

## Author

Sebastian Gregoricchio
