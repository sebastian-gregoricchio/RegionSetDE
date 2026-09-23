# .splitConsensus

Assigns every region of the total consensus to the first set of the user
it overlaps, the others to the unassigned set.

## Usage

``` r
.splitConsensus(
  totalConsensus,
  userList,
  unassignedSet = "other",
  verbose = TRUE
)
```

## Arguments

- totalConsensus:

  `GRanges` with the total consensus.

- userList:

  Named list of `GRanges` with the regions of the user.

- unassignedSet:

  String with the name of the set collecting the regions overlapping no
  set, or `NULL` to drop them.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A named list of `GRanges`, the empty sets left out.

## Author

Sebastian Gregoricchio
