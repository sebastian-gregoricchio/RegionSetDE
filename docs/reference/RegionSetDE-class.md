# RegionSetDE class

S4 class collecting a group of genomic region sets together with the
filters applied to them. The regions are stored as a `GRangesList`, so
that any Bioconductor operation remains available through the `regions`
slot.

## Slots

- `regions`:

  `GRangesList` containing the region sets, one element per set.

- `consensus`:

  List with the consensus data when the regions were built from peaks by
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md):
  the consensus of every group, the peaks of every sample after the
  exclusion, the total consensus, the table of the samples and the
  sample sheet. Empty for regions loaded any other way. Read it with
  [`consensusData`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/consensusData.md).

## Author

Sebastian Gregoricchio
