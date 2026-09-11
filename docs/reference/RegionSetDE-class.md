# RegionSetDE class

S4 class collecting a group of genomic region sets together with the
filters applied to them. The regions are stored as a `GRangesList`, so
that any Bioconductor operation remains available through the `regions`
slot.

## Slots

- `regions`:

  `GRangesList` containing the region sets, one element per set.

## Author

Sebastian Gregoricchio
