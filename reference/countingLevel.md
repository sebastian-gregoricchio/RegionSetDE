# countingLevel

Tells whether the rows of an object are whole regions or tiles of a
region.

## Usage

``` r
countingLevel(object)

# S4 method for class 'RegionSetDE.counts'
countingLevel(object)

# S4 method for class 'RegionSetDE.fit'
countingLevel(object)

# S4 method for class 'RegionSetDE.results'
countingLevel(object)
```

## Arguments

- object:

  `RegionSetDE.counts`, `RegionSetDE.fit` or `RegionSetDE.results`
  object.

## Value

Either `"region"` or `"tile"`.

## See also

[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
[`tileTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/tileTable.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
countingLevel(counts)
#> [1] "region"

fit <- loadExampleData("fit", verbose = FALSE)
countingLevel(fit)
#> [1] "region"
```
