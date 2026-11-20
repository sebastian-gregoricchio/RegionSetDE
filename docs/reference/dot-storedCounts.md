# .storedCounts

Internal accessor of the counts carried by a fit or a result. Unlike
[`resultCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultCounts.md)
it does not stop when the result was built without them, so that the
functions of the package can check for them first.

## Usage

``` r
.storedCounts(object)
```

## Arguments

- object:

  `RegionSetDE.fit`, `RegionSetDE.results` or `RegionSetDE.setResults`
  object.

## Value

A `RegionSetDE.counts` object, with no row and no column when no counts
were carried.

## Author

Sebastian Gregoricchio
