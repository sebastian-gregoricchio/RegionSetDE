# .warnWithinSetOverlap

Raises one warning listing every set with regions overlapping each
other, when the policy asks for it.

## Usage

``` r
.warnWithinSetOverlap(setTable, setColumns, countColumns, policy)
```

## Arguments

- setTable:

  Data.frame of the results, one row per set or per pair.

- setColumns:

  Character vector with the columns naming the sets.

- countColumns:

  Character vector with the columns holding the counts, in the order of
  `setColumns`.

- policy:

  String with the value of `overlapWithinSet`.

## Value

Nothing, called for the warning.

## Author

Sebastian Gregoricchio
