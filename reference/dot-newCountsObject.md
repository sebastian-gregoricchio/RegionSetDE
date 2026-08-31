# .newCountsObject

Assembles a `RegionSetDE.counts` object from a matrix of values, the
regions and the sample table, carrying over the provenance of the region
sets.

## Usage

``` r
.newCountsObject(
  countMatrix,
  regions,
  sampleTable,
  provenance,
  countingLevel = "region",
  newParameters = list(),
  metadataList = list()
)
```

## Arguments

- countMatrix:

  Matrix with the values, one row per region and one column per sample.

- regions:

  `GRanges` used to compute the counts, in the same order as the matrix
  rows.

- sampleTable:

  Data.frame with the sample annotation, in the same order as the matrix
  columns.

- provenance:

  List returned by `.provenanceSlots`.

- countingLevel:

  String indicating whether the rows are regions or sets. Default:
  `"region"`.

- newParameters:

  List with the arguments of the calling function, appended to the
  stored parameters.

- metadataList:

  List stored in the `metadata` of the object. Default:
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A `RegionSetDE.counts` object.

## Author

Sebastian Gregoricchio
