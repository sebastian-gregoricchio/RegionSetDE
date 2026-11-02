# .countsOfObject

Returns the counts object carried by any of the classes holding one, so
that the annotation reads the rows wherever it is called from.

## Usage

``` r
.countsOfObject(object)
```

## Arguments

- object:

  `RegionSetDE.counts`, `RegionSetDE.fit`, `RegionSetDE.results` or
  `RegionSetDE.setResults` object. Not a `RegionSetDE` object: the links
  live beside the counted rows, so the regions have to be counted first.

## Value

A `RegionSetDE.counts` object.

## Author

Sebastian Gregoricchio
