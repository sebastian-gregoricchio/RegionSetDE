# .peakOccupancyClass

Gives every region the class it is summarised under, either the
combination of consensus groups covering it or the number of samples
that carried a peak on it.

## Usage

``` r
.peakOccupancyClass(
  resultTable,
  groups = NULL,
  by = "group",
  contrastGroups = list()
)
```

## Arguments

- resultTable:

  Data.frame with the results, carrying the `peak.*` columns.

- groups:

  Character vector with the consensus groups used, or `NULL`.

- by:

  String, either `"group"` or `"samples"`.

- contrastGroups:

  List with the `column` and the two `groups` of the contrast, possibly
  empty.

## Value

A factor with one value per row of the table.

## Author

Sebastian Gregoricchio
