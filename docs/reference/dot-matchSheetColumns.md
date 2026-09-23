# .matchSheetColumns

Finds the columns of a sample sheet holding the standard fields, from
the names given by the user first, then from the standard names and from
the DiffBind ones, ignoring the case.

## Usage

``` r
.matchSheetColumns(columnNames, columns = NULL)
```

## Arguments

- columnNames:

  Character vector with the column names of the table.

- columns:

  Named character vector given by the user, or `NULL`.

## Value

A data.frame with the `field` and the `column` holding it, one row per
field found.

## Author

Sebastian Gregoricchio
