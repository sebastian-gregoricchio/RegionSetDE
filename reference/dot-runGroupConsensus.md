# .runGroupConsensus

Builds the consensus of one group through
[`consensusRegions::runConsensus`](https://rdrr.io/pkg/consensusRegions/man/runConsensus.html),
called by the workers of
[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

## Usage

``` r
.runGroupConsensus(groupArguments)
```

## Arguments

- groupArguments:

  List with the arguments of
  [`consensusRegions::runConsensus`](https://rdrr.io/pkg/consensusRegions/man/runConsensus.html)
  for the group, the peaks included.

## Value

A `consensusRegions` object.

## Author

Sebastian Gregoricchio
