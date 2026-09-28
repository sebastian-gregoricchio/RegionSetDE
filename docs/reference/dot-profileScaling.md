# .profileScaling

Returns the divisor turning the coverage of every sample into coverage
per million fragments of the average library, scaled by the
normalisation of the object when there is one.

## Usage

``` r
.profileScaling(counts, useOffsets, verbose)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- useOffsets:

  Logical value indicating whether the scaling factors of the
  normalisation must be used.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A numeric vector with one divisor per sample.

## Author

Sebastian Gregoricchio
