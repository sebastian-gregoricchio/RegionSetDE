# .binnedSignal

Reads the signal of every file over a set of windows and averages it in
bins of fixed width, one matrix per file.

## Usage

``` r
.binnedSignal(files, type, windows, binWidth, counts, nThreads = 1)
```

## Arguments

- files:

  Character vector with the BAM or bigWig files, one per sample.

- type:

  String, either `"bam"` or `"bigwig"`.

- windows:

  `GRanges` with the windows, all of the same width, a multiple of
  `binWidth`.

- binWidth:

  Integer value with the width of the bins.

- counts:

  `RegionSetDE.counts` object, whose counting parameters are used for
  BAM files.

- nThreads:

  Number of threads, one file per thread.

## Value

A list with `matrices`, one matrix per file with one row per window and
one column per bin, and `kept`, a logical vector telling which windows
lie within their chromosome.

## Author

Sebastian Gregoricchio
