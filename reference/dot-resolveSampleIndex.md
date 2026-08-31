# .resolveSampleIndex

Turns a sample name or position into a column index.

## Usage

``` r
.resolveSampleIndex(sample = NULL, sampleNames)
```

## Arguments

- sample:

  String or numeric position identifying a sample. Default: `NULL`.

- sampleNames:

  Character vector with the sample names, in the order of the columns.

## Value

The column index, or `NULL` when no sample is given.

## Author

Sebastian Gregoricchio
