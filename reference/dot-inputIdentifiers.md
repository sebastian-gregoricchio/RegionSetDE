# .inputIdentifiers

Names every distinct input library of a sample sheet, checking that
names and files correspond one to one when the names come with the
table.

## Usage

``` r
.inputIdentifiers(inputPaths, inputIds = NULL)
```

## Arguments

- inputPaths:

  Character vector with the input file of every sample, `NA` for the
  samples without one.

- inputIds:

  Character vector with the names given in the table, or `NULL`.

## Value

A character vector with the name of the input of every sample, `NA`
where there is no input.

## Author

Sebastian Gregoricchio
