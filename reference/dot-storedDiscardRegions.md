# .storedDiscardRegions

Returns the regions whose reads `countReads` discarded, written in the
chromosome names of the files about to be read, so that every later
reading of the same BAM files ignores the same reads.

## Usage

``` r
.storedDiscardRegions(counts, targetSeqlevels)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- targetSeqlevels:

  Character vector with the chromosome names of the files.

## Value

A `GRanges`, or `NULL` when no read was discarded at the counting step.

## Author

Sebastian Gregoricchio
