# .placeBrackets

Gives each bracket its horizontal span on the axis and a height above
the data. A bracket shares a level with the others unless the two would
overlap, in which case it moves up one level.

## Usage

``` r
.placeBrackets(bracketTable, groupLevels, valueRange, labelLines = 1)
```

## Arguments

- bracketTable:

  Data.frame with the `group1` and `group2` columns.

- groupLevels:

  Character vector with the groups, in the order of the axis.

- valueRange:

  Numeric vector with the lowest and the highest value drawn.

- labelLines:

  Numeric value with the number of lines of each label. Default: `1`.

## Value

The input data.frame with the `x.start`, `x.end`, `y.bracket`, `y.tick`
and `y.top` columns added.

## Author

Sebastian Gregoricchio
