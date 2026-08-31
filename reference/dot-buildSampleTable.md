# .buildSampleTable

Assembles the sample table used as `colData`, deriving the sample names
from the file names when they are not provided and attaching the user
metadata.

## Usage

``` r
.buildSampleTable(
  files,
  sampleNames = NULL,
  sampleMetadata = NULL,
  fileColumn = "file",
  extensionPattern = "\\.[^.]*$"
)
```

## Arguments

- files:

  Character vector with the paths of the signal files.

- sampleNames:

  Character vector with the sample names. Default: `NULL`, derived from
  the file names.

- sampleMetadata:

  Data.frame with the sample annotation. Default: `NULL`.

- fileColumn:

  String with the name of the column storing the file paths. Default:
  `"file"`.

- extensionPattern:

  Regular expression removed from the file names to build the sample
  names. Default: `"\.[^.]*$"`.

## Value

A data.frame with one row per sample.

## Author

Sebastian Gregoricchio
