# .extraRowColumns

Picks the annotation columns of the `rowData` that must travel into a
result, leaving out the ones the package writes itself and renaming any
that would collide with a statistic.

## Usage

``` r
.extraRowColumns(
  rowTable,
  extraColumns = TRUE,
  reserved = character(0),
  verbose = TRUE
)
```

## Arguments

- rowTable:

  Data.frame with the `rowData` of the counts.

- extraColumns:

  `TRUE` for every annotation column, `FALSE` for none, or a character
  vector naming the ones wanted.

- reserved:

  Character vector with the names already spoken for.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A data.frame with one row per row of the counts, possibly with no column
at all.

## Author

Sebastian Gregoricchio
