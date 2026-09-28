# .consensusWorkers

Shares the threads between the groups built in parallel and the
calibration run inside each group.

## Usage

``` r
.consensusWorkers(nThreads, nGroups, calibrate = FALSE, userBPPARAM = NULL)
```

## Arguments

- nThreads:

  Number of threads asked for.

- nGroups:

  Number of groups needing a consensus, the groups of a single sample
  left out.

- calibrate:

  Logical value to indicate whether the threshold is calibrated inside
  each group.

- userBPPARAM:

  The `BPPARAM` passed by the user to
  [`consensusRegions::runConsensus`](https://rdrr.io/pkg/consensusRegions/man/runConsensus.html),
  or `NULL`.

## Value

A list with `outer`, the number of groups run at once, and `inner`, the
number of threads each group gets for its permutations.

## Author

Sebastian Gregoricchio
