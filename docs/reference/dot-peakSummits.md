# .peakSummits

Places the summit of every region from the peak calls overlapping it:
the summits written in the narrowPeak files, averaged with weights
proportional to their significance.

## Usage

``` r
.peakSummits(regions, peakList)
```

## Arguments

- regions:

  `GRanges` with the regions.

- peakList:

  `GRangesList` with the peaks of every sample, as kept by
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

## Value

An integer vector with the summit of every region, `NA` where no peak
overlaps it.

## Author

Sebastian Gregoricchio
