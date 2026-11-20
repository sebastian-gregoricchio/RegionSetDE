# .joinTechnicalReplicates

Joins the libraries of every biological sample into one peak set through
[`consensusRegions::runConsensus`](https://rdrr.io/pkg/consensusRegions/man/runConsensus.html),
before the consensus of the groups of
[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

## Usage

``` r
.joinTechnicalReplicates(
  sampleTable,
  peakList,
  technicalArguments,
  seqlevelsStyle,
  BPPARAM = NULL,
  verbose = TRUE
)
```

## Arguments

- sampleTable:

  Data.frame with one row per library and the `sample`, `group` and
  `biological` columns.

- peakList:

  Named `GRangesList` with the peaks of every library, after the
  blacklist and the greylist.

- technicalArguments:

  List with the arguments of
  [`consensusRegions::runConsensus`](https://rdrr.io/pkg/consensusRegions/man/runConsensus.html)
  for this step.

- seqlevelsStyle:

  String with the chromosome naming style of the peaks, or `NULL`.

- BPPARAM:

  A `BiocParallelParam` object, a number of cores, or `NULL`. Default:
  `NULL`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A list with `peaks`, a named `GRangesList` with one peak set per
biological sample, and `objects`, a named list with the
`consensusRegions` object behind each of them, `NULL` for the biological
samples with a single library.

## Author

Sebastian Gregoricchio
