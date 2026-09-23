# genomeAssembly

Returns the genome assembly declared for the regions of a RegionSetDE
object.

## Usage

``` r
genomeAssembly(object)

# S4 method for class 'RegionSetDE.provenance'
genomeAssembly(object)
```

## Arguments

- object:

  Any object of the package: `RegionSetDE`, `RegionSetDE.counts`,
  `RegionSetDE.fit`, `RegionSetDE.results`, `RegionSetDE.setResults` or
  `RegionSetDE.setScores`.

## Value

A string with the assembly, or `NULL` when none was declared.

## See also

[`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
genomeAssembly(counts)
#> [1] "rn4"
```
