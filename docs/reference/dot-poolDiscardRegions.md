# .poolDiscardRegions

Pools the regions whose reads must be ignored into one set of ranges
without overlaps, written in the chromosome names of the BAM files. The
lists may come in different naming styles, each of them is brought to
the files on its own.

## Usage

``` r
.poolDiscardRegions(regionLists, targetSeqlevels)
```

## Arguments

- regionLists:

  List of `GRanges`, `NULL` for the lists that were not given.

- targetSeqlevels:

  Character vector with the chromosome names of the BAM files.

## Value

A `GRanges` without overlaps, or `NULL` when no list holds a region.

## Author

Sebastian Gregoricchio
