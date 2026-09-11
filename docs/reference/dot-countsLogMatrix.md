# .countsLogMatrix

Turns a counts object into log2 counts per million, with or without the
normalisation it carries.

## Usage

``` r
.countsLogMatrix(counts, useOffsets = TRUE, priorCount = 2)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- useOffsets:

  Logical value to indicate whether the stored normalisation must be
  applied.

- priorCount:

  Numeric value with the prior count added before taking the logarithm.
  Default: `2`.

## Value

A numeric matrix with one row per region and one column per sample.

## Author

Sebastian Gregoricchio
