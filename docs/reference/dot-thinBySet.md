# .thinBySet

Thins a table of regions down to a number of points that a panel can
hold, inside each region set when the plot is faceted.

## Usage

``` r
.thinBySet(regionTable, maxPoints = 20000, bySet = TRUE)
```

## Arguments

- regionTable:

  Data.frame with a `region.set` column.

- maxPoints:

  Numeric value with the number of points kept per panel.

- bySet:

  Logical value indicating whether the thinning must happen inside each
  set.

## Value

A data.frame with at most `maxPoints` rows per panel.

## Author

Sebastian Gregoricchio
