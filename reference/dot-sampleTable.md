# .sampleTable

Builds the table of the points of an ordination, with the groups driving
the colour and the shape and the labels.

## Usage

``` r
.sampleTable(
  sampleAnnotation,
  colourBy = NULL,
  shapeBy = NULL,
  labelBy = "sample"
)
```

## Arguments

- sampleAnnotation:

  Data.frame with a `sample` column and the annotation of the samples.

- colourBy:

  String with a column of the annotation, or `NULL`.

- shapeBy:

  String with a column of the annotation, or `NULL`.

- labelBy:

  String with a column of the annotation, `"sample"`, or `NULL`.

## Value

A data.frame with the `sample`, `colour.group`, `shape.group` and
`point.label` columns.

## Author

Sebastian Gregoricchio
