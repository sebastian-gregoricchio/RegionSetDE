# .regionAbundance

Computes the average abundance of every row, optionally brought to a
common width so that regions of different sizes can be compared to the
same threshold.

## Usage

``` r
.regionAbundance(
  countMatrix,
  librarySizes,
  rowWidths = NULL,
  referenceWidth = NULL,
  priorCount = 2
)
```

## Arguments

- countMatrix:

  Numeric matrix of counts.

- librarySizes:

  Numeric vector with the library sizes.

- rowWidths:

  Numeric vector with the width of every row, or `NULL` to leave the
  abundances on their own width.

- referenceWidth:

  Numeric value with the width the abundances are scaled to.

- priorCount:

  Numeric value with the prior count added before taking the logarithm.
  Default: `2`.

## Value

A numeric vector of log2 CPM values.

## Author

Sebastian Gregoricchio
