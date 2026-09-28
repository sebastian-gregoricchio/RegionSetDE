# consensusGroupList

Returns the consensus of every group of samples built by
[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md),
as a named list with one `GRanges` per group, ready to be annotated or
compared group by group with tools taking a list of peak sets, such as
`ChIPseeker`.

## Usage

``` r
consensusGroupList(
  object,
  groups = NULL,
  seqlevelsStyle = NULL,
  asGRangesList = FALSE
)
```

## Arguments

- object:

  `RegionSetDE` object returned by
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

- groups:

  Character vector with the groups returned, in the order wanted.
  Default: `NULL`, every group, in the order of the consensus.

- seqlevelsStyle:

  String with the chromosome naming style of the output, one among
  `"UCSC"`, `"Ensembl"` and `"NCBI"`, for instance `"UCSC"` to match a
  `TxDb` of the UCSC annotation. Default: `NULL`, the style of the
  object.

- asGRangesList:

  Logical value to indicate whether a `GRangesList` must be returned
  instead of a list. Default: `FALSE`.

## Value

A named list of `GRanges`, or a `GRangesList`, with one element per
group holding its consensus regions, named after the group.

## Details

The consensus of a group is the one built within that group, before the
groups are pooled, so a region found in two groups is in both elements
and a region found in one group only is in that one alone. The regions
carry no metadata column, since the statistics of the consensus of one
group mean nothing for the others; the consensus object of every group,
with its statistics, is in `consensusData(object)$objects`. A group with
a single sample holds the peaks of that sample, merged where they
overlap, and the peaks removed by the `blacklist` and the `greylist` are
absent from every group.

The list goes as it is to the functions of `ChIPseeker` that take
several peak sets, for instance
`lapply(consensusGroupList(x, seqlevelsStyle = "UCSC"), ChIPseeker::annotatePeak, TxDb = txdb)`
followed by `ChIPseeker::plotAnnoBar()`. The naming style of the
chromosomes has to match the one of the annotation, which is what
`seqlevelsStyle` is for.

## See also

[`consensusData`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/consensusData.md),
[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md),
[`plotPeakUpset`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotPeakUpset.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
if (requireNamespace("consensusRegions", quietly = TRUE)) {
  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
  consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl", verbose = FALSE)

  groupList <- consensusGroupList(consensus, seqlevelsStyle = "UCSC")
  lengths(groupList)
  groupList$R1881_24h

  # Two groups only, in the order given
  names(consensusGroupList(consensus, groups = c("R1881_24h", "DMSO")))
}
#> [1] "R1881_24h" "DMSO"     
```
