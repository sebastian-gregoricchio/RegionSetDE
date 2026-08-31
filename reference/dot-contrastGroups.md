# .contrastGroups

Works out which variable of the sample metadata a contrast separates,
and which two of its levels, by comparing the contrast against the
difference between the design rows of every pair of levels.

## Usage

``` r
.contrastGroups(contrastVector, design, colData = NULL, maxLevels = 20)
```

## Arguments

- contrastVector:

  Numeric vector with the contrast.

- design:

  Design matrix.

- colData:

  Sample metadata, or `NULL`.

- maxLevels:

  Numeric value with the number of levels above which a column is not
  considered a grouping variable. Default: `20`.

## Value

A list with the `column` and the two `groups`, the first one being the
level the contrast is positive for. An empty list when the contrast is
not a difference between two levels of one variable.

## Author

Sebastian Gregoricchio
