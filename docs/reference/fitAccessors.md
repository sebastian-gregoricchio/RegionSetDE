# Fit accessors

Return the parts of a `RegionSetDE.fit` object: what the engine
returned, the engine itself, the design and its formula, the blocking,
the samples, and the universe built for the set level tests. `fitEngine`
and `setUniverse` answer for the results as well, since both travel from
the fit into them.

## Usage

``` r
fitComponents(object)

# S4 method for class 'RegionSetDE.fit'
fitComponents(object)

fitEngine(object)

# S4 method for class 'RegionSetDE.fit'
fitEngine(object)

# S4 method for class 'RegionSetDE.results'
fitEngine(object)

# S4 method for class 'RegionSetDE.setResults'
fitEngine(object)

# S4 method for class 'RegionSetDE.geneResults'
fitEngine(object)

designMatrix(object)

# S4 method for class 'RegionSetDE.fit'
designMatrix(object)

designFormula(object)

# S4 method for class 'RegionSetDE.fit'
designFormula(object)

blockingInfo(object)

# S4 method for class 'RegionSetDE.fit'
blockingInfo(object)

fitSamples(object)

# S4 method for class 'RegionSetDE.fit'
fitSamples(object)

setUniverse(object)

# S4 method for class 'RegionSetDE.fit'
setUniverse(object)

# S4 method for class 'RegionSetDE.setResults'
setUniverse(object)
```

## Arguments

- object:

  `RegionSetDE.fit` object returned by
  [`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md).
  `fitEngine` also takes a `RegionSetDE.results`,
  `RegionSetDE.setResults` or `RegionSetDE.geneResults` object, and
  `setUniverse` a `RegionSetDE.setResults` object.

## Value

`fitComponents` returns a named list with everything the engine
produced: `object`, the fit itself, which
[`fitObject`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitObject.md)
returns alone, and the auxiliary objects needed to test it, such as the
voom weights or the sample table. `fitEngine` returns a string, one
among `"edgeR"`, `"voom"`, `"dream"`, `"deseq2"` and `"limma"`.
`designMatrix` returns the matrix of the fixed effects, one row per
sample, as `design` does. `designFormula` returns the formula of the
design, or `NULL` when the design was given as a matrix. `blockingInfo`
returns a named list with the blocking variable, the random effect
formula and the consensus correlation, empty when there is no blocking.
`fitSamples` returns the names of the samples the model was fitted on.
`setUniverse` returns the `RegionSetDE.universe` object with the rows
every set is compared against, read with the [universe
accessors](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/universeAccessors.md).

## See also

[`fitObject`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitObject.md),
[`fitCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitCounts.md),
[`dispersionInfo`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/dispersionInfo.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

fitEngine(fit)
#> [1] "edgeR"
designFormula(fit)
#> ~condition
#> <environment: base>
head(designMatrix(fit))
#>                                 (Intercept) conditionSHR
#> lv-H3K4me3-BN-female-bio1-tech1           1            0
#> lv-H3K4me3-BN-male-bio2-tech1             1            0
#> lv-H3K4me3-SHR-male-bio2-tech1            1            1
#> lv-H3K4me3-SHR-male-bio3-tech1            1            1
fitSamples(fit)
#> [1] "lv-H3K4me3-BN-female-bio1-tech1" "lv-H3K4me3-BN-male-bio2-tech1"  
#> [3] "lv-H3K4me3-SHR-male-bio2-tech1"  "lv-H3K4me3-SHR-male-bio3-tech1" 
names(fitComponents(fit))
#> [1] "object" "dge"    "test"  
blockingInfo(fit)
#> list()
setUniverse(fit)
#> An object of class 'RegionSetDE.universe'
#>   type            : otherSets 
#>   matched on      : width, abundance 
#>   compared against: promoterNonCpG, intergenic, geneBody, promoterCpG 
#>   sets            : 4 
#> 
#>      region.set n.regions n.comparison median.width median.width.comparison
#>  promoterNonCpG       277         1378         1000                    1000
#>      intergenic       440         1354         1000                    1000
#>        geneBody       909          986         1000                    1000
#>     promoterCpG       269          314         1000                    1000
#>  median.abundance median.abundance.comparison
#>              3.77                        3.78
#>              3.44                        3.76
#>              3.60                        4.08
#>              9.81                        3.75
```
