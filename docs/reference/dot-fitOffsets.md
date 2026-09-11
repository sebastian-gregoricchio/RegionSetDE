# .fitOffsets

Extracts the normalisation stored in a counts object as a matrix of log
offsets.

## Usage

``` r
.fitOffsets(counts, useOffsets = TRUE, verbose = TRUE)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- useOffsets:

  Logical value indicating whether the normalisation must be used.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A numeric matrix of log offsets, or `NULL`.

## Author

Sebastian Gregoricchio
