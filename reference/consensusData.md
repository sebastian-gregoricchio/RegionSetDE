# consensusData

Returns the consensus data kept by a `RegionSetDE` object built from
peaks with
[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

## Usage

``` r
consensusData(object)

# S4 method for class 'RegionSetDE'
consensusData(object)
```

## Arguments

- object:

  `RegionSetDE` object.

## Value

A list with `groups`, the consensus regions of every group; `objects`,
the `consensusRegions` object behind each of them (`NULL` for the groups
with a single sample); `total`, the pooled consensus; `peaks`, the peaks
of every sample after the blacklist and the greylist; `removed`, the
peaks they took out, with the `sample` they came from and the list that
removed them in `removed.by`; `samples`, a data.frame with the group,
the peak file, the number of peaks read (`n.peaks`), removed by the
blacklist (`n.blacklist`), by the greylist (`n.greylist`) and by both
(`n.excluded`) for every sample; `sheet`, the sample sheet the consensus
was built from; `groupBy`; `mode`, one among `"consensus"`, `"split"`
and `"replace"`; and `blacklist` and `greylist`, the regions the peaks
were cleaned against, `NULL` when a list was not given.

## See also

[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md),
[`consensusGroupList`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/consensusGroupList.md),
[`plotPeakUpset`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotPeakUpset.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
if (requireNamespace("consensusRegions", quietly = TRUE)) {
  peakFiles <- list.files(system.file("extdata", package = "consensusRegions"),
                          pattern = "rep[0-9]\\.narrowPeak$", full.names = TRUE)

  sampleSheet <- loadSampleSheet(data.frame(sample = c("A_1", "A_2", "B_1", "B_2"),
                                            bam = c("A_1.bam", "A_2.bam", "B_1.bam", "B_2.bam"),
                                            peaks = peakFiles[c(1, 2, 2, 3)],
                                            condition = c("A", "A", "B", "B")),
                                 checkFiles = FALSE, verbose = FALSE)

  regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)

  consensusData(regions)$samples
  lengths(consensusData(regions)$groups)
}
#>   A   B 
#> 287 262 
```
