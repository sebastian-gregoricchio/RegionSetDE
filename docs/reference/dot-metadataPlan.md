# .metadataPlan

Works out which metadata columns of a collection of region sets can be
stacked into one table, and in which type, so that sets loaded from
different files can be put one after the other.

## Usage

``` r
.metadataPlan(regionList, keepMetadata = TRUE, verbose = TRUE)
```

## Arguments

- regionList:

  Named list of `GRanges`.

- keepMetadata:

  Logical value to indicate whether the metadata must be carried over at
  all.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A list with, for every column kept, the name it takes in the output and
the type it is stored as.

## Author

Sebastian Gregoricchio
