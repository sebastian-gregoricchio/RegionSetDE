# .asBiologicalPeaks

Brings the peaks of a biological sample to the columns the consensus of
a group reads: `pValue`, the significance on the -log10 scale, and
`peak`, the summit as an offset from the start, as in a narrowPeak file.

## Usage

``` r
.asBiologicalPeaks(peakRanges, libraryPeaks = NULL)
```

## Arguments

- peakRanges:

  `GRanges` with the peaks of a single library, carrying `negLog10P`, or
  with the consensus of several libraries, carrying `combinedNegLog10P`.

- libraryPeaks:

  Named `GRangesList` with the peaks of the libraries behind a
  consensus, used to place its summits, or `NULL` for a single library.

## Value

A `GRanges` with the `pValue` and `peak` columns only. `peak` is -1
where the peaks carry no summit.

## Author

Sebastian Gregoricchio
