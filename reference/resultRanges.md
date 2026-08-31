# resultRanges

Returns the coordinates of a `RegionSetDE.results` object with the
statistics attached as metadata columns, ready to be exported as a
BED-like file.

## Usage

``` r
resultRanges(results)

# S4 method for class 'RegionSetDE.results'
resultRanges(results)

# S4 method for class 'RegionSetDE.resultsList'
resultRanges(results)
```

## Arguments

- results:

  `RegionSetDE.results` object.

## Value

A `GRanges` with one element per region.

## Author

Sebastian Gregoricchio
