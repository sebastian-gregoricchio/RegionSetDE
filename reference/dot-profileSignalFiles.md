# .profileSignalFiles

Works out which files the signal of a profile is read from, and of which
type they are.

## Usage

``` r
.profileSignalFiles(counts, signal, signalFiles)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- signal:

  String, one among `"auto"`, `"bam"` and `"bigwig"`.

- signalFiles:

  Character vector with one file per sample, or `NULL`.

## Value

A list with the `files` and their `type`, either `"bam"` or `"bigwig"`.

## Author

Sebastian Gregoricchio
