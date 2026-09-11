# .pickResults

Returns the single `RegionSetDE.results` object a function has to work
on, whether it was handed one directly or one contrast out of a list.

## Usage

``` r
.pickResults(results, contrast = NULL)
```

## Arguments

- results:

  `RegionSetDE.results`, `RegionSetDE.setResults`, or one of the two
  list classes holding several of them.

- contrast:

  String with the name of the contrast, or its position. Default:
  `NULL`.

## Value

A `RegionSetDE.results` or a `RegionSetDE.setResults` object.

## Author

Sebastian Gregoricchio
