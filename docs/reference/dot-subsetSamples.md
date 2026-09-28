# .subsetSamples

Keeps the samples of a counts object selected for a figure, together
with the normalisation stored for them.

## Usage

``` r
.subsetSamples(counts, samples = NULL)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- samples:

  Character, numeric or logical vector with the samples kept, or `NULL`
  for all of them.

## Value

The `RegionSetDE.counts` object restricted to the samples kept.

## Author

Sebastian Gregoricchio
