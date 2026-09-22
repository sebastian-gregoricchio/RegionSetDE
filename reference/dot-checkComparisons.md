# .checkComparisons

Checks the pairs of groups asked for through `comparisons`, or builds
every pair when none is given.

## Usage

``` r
.checkComparisons(comparisons = NULL, groupLevels)
```

## Arguments

- comparisons:

  List of character vectors of length two, or `NULL`. Default: `NULL`.

- groupLevels:

  Character vector with the groups, in the order of the axis.

## Value

A list of character vectors of length two, each pair appearing once.

## Author

Sebastian Gregoricchio
