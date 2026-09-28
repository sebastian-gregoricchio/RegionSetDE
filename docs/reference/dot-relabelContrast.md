# .relabelContrast

Applies new cut-offs to the result of a single contrast: new labels for
the regions, new thresholds in the object and in its parameters.

## Usage

``` r
.relabelContrast(results, FDR, log2FC, verbose)
```

## Arguments

- results:

  `RegionSetDE.results` or `RegionSetDE.setResults` object.

- FDR:

  Numeric value with the new adjusted p-value cut-off, or `NULL` to keep
  the stored one.

- log2FC:

  Numeric value with the new absolute log2 fold change cut-off, or
  `NULL` to keep the stored one.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

The object with the new labels and thresholds.

## Author

Sebastian Gregoricchio
