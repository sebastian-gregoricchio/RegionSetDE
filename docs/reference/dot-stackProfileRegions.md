# .stackProfileRegions

Stacks the groups of rows of a profile into one `GRanges`, dropping
their seqinfo so that groups from different sources can be pasted
together.

## Usage

``` r
.stackProfileRegions(groupList)
```

## Arguments

- groupList:

  List of `GRanges`, one per group of rows.

## Value

A `GRanges`.

## Author

Sebastian Gregoricchio
