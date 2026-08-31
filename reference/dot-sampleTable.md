# .sampleTable

Assembles the sample metadata a sample level plot maps onto colour,
shape and labels.

## Usage

``` r
.sampleTable(counts, colourBy = NULL, shapeBy = NULL, labelBy = "sample")
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- colourBy:

  String with a `colData` column, or `NULL`.

- shapeBy:

  String with a `colData` column, or `NULL`.

- labelBy:

  String with a `colData` column, `"sample"`, or `NULL`.

## Value

A data.frame with one row per sample.

## Author

Sebastian Gregoricchio
