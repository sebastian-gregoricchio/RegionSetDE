# RegionSetDE.counts class

S4 class storing the read counts computed over a collection of region
sets. It extends `RangedSummarizedExperiment`, therefore `assay`,
`colData`, `rowRanges` and the subsetting operators behave as usual,
while the filters applied upstream remain accessible in the inherited
provenance slots.

## Slots

- `counting.level`:

  String indicating whether the rows are regions (`"region"`) or tiles
  of a region (`"tile"`). It is the same vocabulary the
  `RegionSetDE.fit` and `RegionSetDE.results` classes use, since the
  value travels from here into both of them.

## Author

Sebastian Gregoricchio
