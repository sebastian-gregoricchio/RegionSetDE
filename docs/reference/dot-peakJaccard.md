# .peakJaccard

Computes the Jaccard index between the peak sets of every pair of
samples, on the regions of the union of the peaks or on the base pairs
they cover.

## Usage

``` r
.peakJaccard(object, set = NULL, jaccardLevel = "region", samples = NULL)
```

## Arguments

- object:

  `RegionSetDE` object built by
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md),
  or a named `GRangesList` or list of `GRanges` with the peaks of every
  sample.

- set:

  Character vector with the region sets whose regions the peaks must
  overlap, or `NULL` for all the peaks.

- jaccardLevel:

  String, either `"region"` or `"basepair"`.

- samples:

  Character, numeric or logical vector with the samples kept, or `NULL`
  for all of them.

## Value

A list shaped as the one of
[`computeSampleCorrelation`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeSampleCorrelation.md):
`correlation`, the matrix of the Jaccard indices, `samples`, the
annotation of the samples, and `parameters`.

## Author

Sebastian Gregoricchio
