# .nullMatrix

Collects the counts and the offsets of the rows a null estimate is
computed on.

## Usage

``` r
.nullMatrix(
  counts,
  source = "background",
  regionSets = NULL,
  index = NULL,
  subset = NULL,
  minCount = 10,
  maxRows = 50000
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- source:

  String with where the null rows come from.

- regionSets:

  Character vector with the sets used as null rows, or `NULL`.

- index:

  Integer vector with the positions of the null rows, or `NULL`.

- subset:

  Integer vector restricting the selected rows further, or `NULL`.

- minCount:

  Numeric value with the average count a row must carry.

- maxRows:

  Numeric value with the number of rows kept.

## Value

A list with the filtered `counts` and `offset` matrices, the
`library.size` vector, the `abundance` of the kept rows and `kept.rows`,
their positions in the unfiltered selection.

## Author

Sebastian Gregoricchio
