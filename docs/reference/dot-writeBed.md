# .writeBed

Writes the coordinates of a result table as a BED file, coloured by the
direction of the change when asked for.

## Usage

``` r
.writeBed(resultTable, path, bedScore = "FDR", colourByStatus = TRUE)
```

## Arguments

- resultTable:

  Data.frame with the `seqnames`, `start` and `end` columns.

- path:

  String with the destination.

- bedScore:

  String with the column mapped onto the score.

- colourByStatus:

  Logical value indicating whether an item colour must be written.

## Value

Invisibly the path.

## Author

Sebastian Gregoricchio
