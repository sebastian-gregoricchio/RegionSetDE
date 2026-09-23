# .annotationColourList

Gives the colours of the annotation columns of a heatmap, from the user
when given, a grey gradient for the numeric columns, and a palette per
categorical column otherwise, the first one being the palette the other
plots of the package give to the groups.

## Usage

``` r
.annotationColourList(sampleTable, annotationColumns, annotationColours = NULL)
```

## Arguments

- sampleTable:

  Data.frame with the annotation of the samples.

- annotationColumns:

  Character vector with the columns drawn.

- annotationColours:

  Named list with the colours given by the user, or `NULL`.

## Value

A named list with the colours of every annotation column, as
`ComplexHeatmap` takes them.

## Author

Sebastian Gregoricchio
