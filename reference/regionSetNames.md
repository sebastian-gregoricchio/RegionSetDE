# regionSetNames

Returns the names of the region sets carried by an object, at any stage
of the analysis.

## Usage

``` r
regionSetNames(object)

# S4 method for class 'RegionSetDE'
regionSetNames(object)

# S4 method for class 'RegionSetDE.counts'
regionSetNames(object)

# S4 method for class 'RegionSetDE.fit'
regionSetNames(object)

# S4 method for class 'RegionSetDE.results'
regionSetNames(object)

# S4 method for class 'RegionSetDE.setResults'
regionSetNames(object)

# S4 method for class 'RegionSetDE.resultsList'
regionSetNames(object)

# S4 method for class 'RegionSetDE.setResultsList'
regionSetNames(object)
```

## Arguments

- object:

  `RegionSetDE`, `RegionSetDE.counts`, `RegionSetDE.fit` or
  `RegionSetDE.results` object.

## Value

Character vector with the names of the region sets, in the order they
are stored.

## See also

[`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md),
[`splitLoadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/splitLoadRegions.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
regionSetNames(counts)
#> [1] "promoterNonCpG" "intergenic"     "geneBody"       "promoterCpG"   
```
