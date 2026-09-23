# .checkSheetFiles

Checks that the files of a sample sheet exist, and that the BAM files
are indexed, reporting every problem at once.

## Usage

``` r
.checkSheetFiles(sheetTable, pathFields)
```

## Arguments

- sheetTable:

  Data.frame with the standard fields already resolved.

- pathFields:

  Character vector with the fields holding paths.

## Value

Nothing, it stops when a file is missing or a BAM file has no index.

## Author

Sebastian Gregoricchio
