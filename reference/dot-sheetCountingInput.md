# .sheetCountingInput

Takes the signal files, the sample names and the annotation of the
counting from a sample sheet, or from the sample sheet a consensus was
built from when no file is given at all.

## Usage

``` r
.sheetCountingInput(
  regionSet,
  sampleSheet,
  files,
  sampleNames,
  sampleMetadata,
  fileField,
  verbose = TRUE
)
```

## Arguments

- regionSet:

  Object handed to the counting function.

- sampleSheet:

  Data.frame returned by `loadSampleSheet`, the path to a sample sheet,
  or `NULL`.

- files:

  Character vector with the files given directly, or `NULL`.

- sampleNames:

  Character vector with the sample names given directly, or `NULL`.

- sampleMetadata:

  Data.frame with the annotation given directly, or `NULL`.

- fileField:

  String with the column of the sheet holding the files, either `"bam"`
  or `"bigwig"`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A list with the `files`, the `sampleNames` and the `sampleMetadata` to
count with.

## Author

Sebastian Gregoricchio
