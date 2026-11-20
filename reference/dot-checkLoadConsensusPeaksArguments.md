# .checkLoadConsensusPeaksArguments

Checks the arguments of
[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md)
and brings them to the form the rest of the function works with.

## Usage

``` r
.checkLoadConsensusPeaksArguments(
  sampleSheet,
  groupBy,
  technicalBy,
  technicalArguments,
  regionMode,
  unassignedSet,
  seed,
  weightMethod,
  calibrate,
  nPermutations,
  targetFDR,
  recentre,
  width,
  combinationMethod,
  replicateType,
  stringencyThreshold,
  weakThreshold,
  combinedThreshold,
  alpha,
  minSupport,
  minReplicates,
  minOverlap,
  minOverlapFraction,
  mergeMethod,
  maxConsensusWidth,
  verbose,
  ...
)
```

## Arguments

- sampleSheet:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- groupBy:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- technicalBy:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- technicalArguments:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- regionMode:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- unassignedSet:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- seed:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- weightMethod:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- calibrate:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- nPermutations:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- targetFDR:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- recentre:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- width:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- combinationMethod:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- replicateType:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- stringencyThreshold:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- weakThreshold:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- combinedThreshold:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- alpha:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- minSupport:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- minReplicates:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- minOverlap:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- minOverlapFraction:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- mergeMethod:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- maxConsensusWidth:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- verbose:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- ...:

  See
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

## Value

A list with the arguments as checked and completed:
`consensusArguments`, `technicalArguments`, the arguments of the step
joining the technical replicates (`NULL` without `technicalBy`),
`regionMode`, `sampleSheet`.

## Author

Sebastian Gregoricchio
