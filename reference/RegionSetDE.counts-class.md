# RegionSetDE.counts class

S4 class storing the read counts computed over a collection of region
sets. It extends `RangedSummarizedExperiment`, therefore `assay`,
`colData`, `rowRanges` and the subsetting operators behave as usual,
while the filters applied upstream remain accessible in the inherited
provenance slots.

## Slots

- `counting.level`:

  String indicating whether the rows correspond to the single regions
  (`"region"`) or to the whole sets (`"set"`).

## Author

Sebastian Gregoricchio
