# loadConsensusPeaks

Builds the regions of an analysis from the peaks called on each sample.
The peaks are combined into a consensus within each group of samples
through `consensusRegions`, and the group consensus are pooled into a
total one, which becomes the region set counted downstream. Regions
supplied by the user can split the total consensus into sets, or take
its place.

## Usage

``` r
loadConsensusPeaks(
  sampleSheet,
  groupBy = NULL,
  blacklist = NULL,
  greylist = NULL,
  regionSets = NULL,
  regionMode = "split",
  unassignedSet = "other",
  seqlevelsStyle = "UCSC",
  genomeAssembly = NULL,
  nThreads = 1,
  seed = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- sampleSheet:

  Data.frame returned by
  [`loadSampleSheet`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadSampleSheet.md),
  or the path to a sample sheet, with at least the `sample` and `peaks`
  columns.

- groupBy:

  String with the column of the sample sheet defining the groups, for
  instance `"condition"`. Default: `NULL`, all the samples form a single
  group.

- blacklist:

  Regions of the assembly whose peaks must be dropped before any
  consensus is built, typically the list returned by
  [`loadBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadBlacklist.md).
  Either a `GRanges`, a path to a BED-like file, a data.frame, or a list
  of them, which are pooled. The list is stored in the `blacklist` slot
  of the object, as
  [`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md)
  does. Default: `NULL`.

- greylist:

  Regions of this experiment whose peaks must be dropped before any
  consensus is built, typically the list built from the inputs by
  [`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md).
  Accepts the same forms as `blacklist`. Default: `NULL`.

- regionSets:

  Regions of interest of the user, in any form accepted by
  [`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md):
  a named list of paths, `GRanges` or data.frames. Default: `NULL`, the
  total consensus is the only set.

- regionMode:

  String indicating what `regionSets` do, either `"split"`, where every
  region of the total consensus goes to the first set it overlaps, or
  `"replace"`, where the sets of the user are the regions and the peaks
  only annotate them. Ignored without `regionSets`. Default: `"split"`.

- unassignedSet:

  String with the name of the set collecting the consensus regions that
  overlap none of `regionSets` in `"split"` mode. `NULL` drops them.
  Default: `"other"`.

- seqlevelsStyle:

  String indicating the chromosome naming style, one among `"UCSC"`,
  `"Ensembl"` or `"NCBI"`, or `NULL` to keep the names as they are.
  Default: `"UCSC"`.

- genomeAssembly:

  String indicating the genome assembly to store with the regions, e.g.
  `"hg38"`. Default: `NULL`.

- nThreads:

  Number of threads. The groups with more than one sample are built side
  by side, one worker per group, and the threads left over go to the
  calibration of the threshold inside each group when `calibrate = TRUE`
  is passed on to
  [`consensusRegions::runConsensus`](https://rdrr.io/pkg/consensusRegions/man/runConsensus.html).
  Default: `1`.

- seed:

  Number seeding the random numbers of the calibration, when
  `calibrate = TRUE` is passed on to
  [`consensusRegions::runConsensus`](https://rdrr.io/pkg/consensusRegions/man/runConsensus.html).
  It goes to `runConsensus` for every group, which hands it to the
  `BiocParallel` back end of the permutations as its `RNGseed`, and it
  is stored in `parameters$loadConsensusPeaks$seed`, so that the
  consensus can be built again identically, on any number of threads.
  Needs consensusRegions 0.99.2 or later. Default: `NULL`, no seed, and
  a calibrated consensus may change from one call to the next.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

- ...:

  Further arguments passed to
  [`consensusRegions::runConsensus`](https://rdrr.io/pkg/consensusRegions/man/runConsensus.html),
  for instance `minReplicates`, `combinedThreshold`, `calibrate` or
  `weightMethod`.

## Value

A `RegionSetDE` object. Every region carries `peak.<group>`, telling
whether the consensus of that group overlaps it, `peak.groups`, the
number of groups with a consensus peak on it, and `peak.samples`, the
number of samples with a peak on it. The `consensus` slot, read with
[`consensusData`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/consensusData.md),
keeps the consensus of every group, the `consensusRegions` object behind
it, the peaks of every sample after the blacklist and the greylist, the
peaks they removed, the two lists, the total consensus, the table of the
samples and the sample sheet itself, which
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
and
[`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md)
read when no file is given to them. The blacklist goes to the
`blacklist` slot, the greylist to the `greylist` slot with its size in
`parameters$greylist`, and the peaks each of them removed from every
sample to the `filtering.log`, read with
[`filteringLog`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filteringLog.md).
Where the random numbers of the calibration came from is recorded in
`parameters$loadConsensusPeaks`, see Details.

## Details

The consensus is built group by group and then pooled, rather than once
over all the samples. A single requirement over all the libraries, such
as a peak in at least two of them, favours the larger group, while the
same requirement applied within each group treats them alike, and the
union keeps the regions found in one group only. A group with a single
sample has no replicate to agree with, and its peaks stand for the group
as they are.

The peaks are read by `consensusRegions`, which keeps the standard
chromosomes only, as
[`GenomeInfoDb::keepStandardChromosomes`](https://rdrr.io/pkg/GenomeInfoDb/man/seqlevels-wrappers.html)
defines them: scaffolds, patches and unplaced contigs are dropped, and
so are chromosomes whose names it does not recognise.

The peaks overlapping `blacklist` and `greylist` are removed before the
consensus, not after it. An artefact lying next to a genuine peak would
otherwise merge with it, and removing the merged region afterwards would
take the genuine peak away as well. The blacklist is applied first and
the greylist to the peaks left, so a peak lying on both is counted as
blacklisted. The two lists are kept apart because they say different
things: a blacklist describes the assembly and is the same for every
experiment, a greylist describes the inputs of this one. Both are
recorded the way
[`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md)
and
[`applyGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyGreylist.md)
record them, so the counts, the fit and the results built from the
object carry the two lists with them, and
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
leaves their reads out of the library sizes. `consensusData(x)$removed`
holds the peaks taken out, with the sample they came from and the list
that removed them, and `consensusData(x)$samples` their number per
sample. A blacklist built for another assembly than `genomeAssembly` is
refused, as
[`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md)
refuses it.

With `regionMode = "split"` the consensus regions are assigned to the
sets of `regionSets` in the order they are given, a region overlapping
two sets going to the first one. The sets then describe classes of
peaks, such as promoter and distal ones, and
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
can compare them. With `regionMode = "replace"` the counted regions are
those of `regionSets`, and the peaks only tell which of them are
occupied in each group.

The occupancy columns describe where the peaks were called, and they are
a poor basis for a set of regions to test. A set of the regions found in
one group only was selected on the signal of that group, so a test of
its change between the groups answers a question already settled by the
selection.

The consensus of a group does not depend on the other groups, so with
`nThreads` above one the groups run in parallel. `consensusRegions`
builds a consensus on a single thread and only parallelises the
permutations of the calibration, so without `calibrate` the groups are
the one place where more threads save time: three groups on three
threads take about as long as the largest of them. With
`calibrate = TRUE` the threads are shared, `nThreads` divided by the
number of groups for the permutations of each group. A `BPPARAM` passed
in `...` is handed to every group as it is and the groups then run one
after the other, so that two levels of workers are never stacked on each
other by accident.

The calibration is the only step drawing random numbers: it places the
peaks at random to see how often they overlap by chance. The draws come
from `BiocParallel`, and without a seed of its own a `BiocParallel` back
end takes them from a stream it keeps for the whole session, which a
[`set.seed()`](https://rdrr.io/r/base/Random.html) before the call does
not reset: the calibrated threshold, and with it the consensus, can then
change from one call to the next. `seed` is passed on to
[`consensusRegions::runConsensus`](https://rdrr.io/pkg/consensusRegions/man/runConsensus.html)
in every group, which gives it to the back end of the permutations as
its `RNGseed`: the same seed returns the same consensus whatever the
number of threads, and the random numbers of the session are left as
they were. Neither package calls
[`set.seed()`](https://rdrr.io/r/base/Random.html).
`parameters$loadConsensusPeaks` keeps the `seed` and its `seed.source`:
`"seed"`, `"BPPARAM"` when a `BPPARAM` passed in `...` carries its own
`RNGseed`, which is then the one recorded, `"unseeded"` when the
calibration ran without a seed, or `"none"` without calibration. A
`BPPARAM` passed in `...` is used as it is, so `seed` is ignored with
it, with a warning.

## See also

[`plotPeakUpset`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotPeakUpset.md),
[`consensusData`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/consensusData.md),
[`loadSampleSheet`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadSampleSheet.md),
[`loadBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadBlacklist.md),
[`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
if (requireNamespace("consensusRegions", quietly = TRUE)) {
  peakFiles <- list.files(system.file("extdata", package = "consensusRegions"),
                          pattern = "rep[0-9]\\.narrowPeak$", full.names = TRUE)

  # Three replicates of one condition and two of another, peaks only
  sampleSheet <- data.frame(sample = c("A_1", "A_2", "A_3", "B_1", "B_2"),
                            bam = c("A_1.bam", "A_2.bam", "A_3.bam", "B_1.bam", "B_2.bam"),
                            peaks = peakFiles[c(1, 2, 3, 1, 2)],
                            condition = c("A", "A", "A", "B", "B"))

  sampleSheet <- loadSampleSheet(sampleSheet, checkFiles = FALSE, verbose = FALSE)

  regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition")
  regions

  head(regionRanges(regions)$consensus, 3)

  # The same peaks split by regions of interest of the user
  firstHalf <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 5e5))
  splitRegions <- loadConsensusPeaks(sampleSheet, groupBy = "condition",
                                     regionSets = list(firstHalf = firstHalf),
                                     verbose = FALSE)
  lengths(regionRanges(splitRegions))
}
#> Group A: 3 samples, 292 consensus regions.
#> Group B: 2 samples, 287 consensus regions.
#> Total consensus: 292 regions from 2 groups.
#> firstHalf     other 
#>         7       285 
```
