# Universe accessors

Return the parts of a `RegionSetDE.universe` object, the rows every
region set is compared against in a competitive test.
[`setUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitAccessors.md)
reads the universe out of a fit or a set level result.

## Usage

``` r
universeIndex(universe)

# S4 method for class 'RegionSetDE.universe'
universeIndex(universe)

universeType(universe)

# S4 method for class 'RegionSetDE.universe'
universeType(universe)

universeMatching(universe)

# S4 method for class 'RegionSetDE.universe'
universeMatching(universe)

comparisonSets(universe)

# S4 method for class 'RegionSetDE.universe'
comparisonSets(universe)

universeDiagnostics(universe)

# S4 method for class 'RegionSetDE.universe'
universeDiagnostics(universe)

universeRows(universe)

# S4 method for class 'RegionSetDE.universe'
universeRows(universe)
```

## Arguments

- universe:

  `RegionSetDE.universe` object, returned by
  [`makeSetUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md)
  or by
  [`setUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitAccessors.md).

## Value

`universeIndex` returns a named list with one element per region set,
the positions of the rows forming its universe. `universeType` returns a
string with the way the universe was built, for instance `"otherSets"`,
`"supplied"`, `"allRows"` or `"pairedSet"`. `universeMatching` returns
the covariates the comparison rows were matched on, empty when they were
not matched. `comparisonSets` returns the names of the sets the
comparison rows came from. `universeDiagnostics` returns a data.frame
with, for every set, the median width and abundance of the set and of
the rows it is compared against. `universeRows` returns the number of
rows of the object the positions refer to.

## See also

[`makeSetUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md),
[`plotUniverseMatching`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotUniverseMatching.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
universe <- setUniverse(fit)

universeType(universe)
#> [1] "otherSets"
universeMatching(universe)
#> [1] "width"     "abundance"
comparisonSets(universe)
#> [1] "promoterNonCpG" "intergenic"     "geneBody"       "promoterCpG"   
lengths(universeIndex(universe))
#> promoterNonCpG     intergenic       geneBody    promoterCpG 
#>           1655           1794           1895            583 
universeDiagnostics(universe)
#>       region.set n.regions n.comparison median.width median.width.comparison
#> 1 promoterNonCpG       277         1378         1000                    1000
#> 2     intergenic       440         1354         1000                    1000
#> 3       geneBody       909          986         1000                    1000
#> 4    promoterCpG       269          314         1000                    1000
#>   median.abundance median.abundance.comparison
#> 1             3.77                        3.78
#> 2             3.44                        3.76
#> 3             3.60                        4.08
#> 4             9.81                        3.75
universeRows(universe)
#> [1] 1895
```
