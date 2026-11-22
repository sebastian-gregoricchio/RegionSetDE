# Set level result accessors

Return the parts of a `RegionSetDE.setResults` object beyond its table:
the per-region statistics the set level test was computed on, the kind
of test and the methods that were run.

## Usage

``` r
regionStats(results)

# S4 method for class 'RegionSetDE.setResults'
regionStats(results)

setTestType(results)

# S4 method for class 'RegionSetDE.setResults'
setTestType(results)

setTestMethods(results)

# S4 method for class 'RegionSetDE.setResults'
setTestMethods(results)
```

## Arguments

- results:

  `RegionSetDE.setResults` object returned by
  [`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
  or
  [`testSetContrast`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md).

## Value

`regionStats` returns a data.frame with one row per region, holding the
statistics the set level test read. `setTestType` returns `"set"` for
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
and `"setContrast"` for
[`testSetContrast`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md).
`setTestMethods` returns the names of the tests that were run, for
instance `"camera"` and `"fry"`.

## See also

[`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md),
[`setUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitAccessors.md),
[`contrastGroups`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/contrastAccessors.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
setRes <- testRegionSets(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

setTestType(setRes)
#> [1] "set"
setTestMethods(setRes)
#> [1] "camera" "fry"   
head(regionStats(setRes), 3)
#>       log2FC average.signal       stat stat.distribution df1      df2   p.value
#> 1 -0.4744876       3.098695 -0.4077872                 f   1 18.94543 0.6880013
#> 2 -1.6911510       3.141973 -1.5289174                 f   1 18.73556 0.1429931
#> 3  0.1328366       3.273341  0.1001436                 f   1 18.45274 0.9213104
#>       region.set    region.id region.key.plain seqnames start   end width
#> 1 promoterNonCpG region_00012     region_00012    chr12 26988 27987  1000
#> 2 promoterNonCpG region_00017     region_00017    chr12 39449 40448  1000
#> 3 promoterNonCpG region_00019     region_00019    chr12 44116 45115  1000
```
