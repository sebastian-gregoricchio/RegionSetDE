# .flattenParameters

Turns the nested list of parameters an object carries into one line per
value, so that the record of an analysis can be read, grepped and
compared against another.

## Usage

``` r
.flattenParameters(x, prefix = "", maxLength = 200)
```

## Arguments

- x:

  List, or any value held inside one.

- prefix:

  String with the path to the current value, built as the recursion
  descends.

- maxLength:

  Numeric value with the number of characters a value is truncated to.
  Default: `200`.

## Value

A data.frame with the `parameter` and `value` columns.

## Author

Sebastian Gregoricchio
