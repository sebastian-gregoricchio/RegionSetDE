# .resolveInputFiles

Returns the input file of every sample, given directly or taken from the
`input` column of the sample annotation.

## Usage

``` r
.resolveInputFiles(inputFiles, sampleTable)
```

## Arguments

- inputFiles:

  Character vector given to `countReads`, or `NULL`.

- sampleTable:

  Data.frame with one row per sample.

## Value

A character vector with one path per sample, `NA` where there is no
input, or `NULL` when no sample has one.

## Author

Sebastian Gregoricchio
