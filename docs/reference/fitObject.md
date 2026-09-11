# fitObject

Returns the object produced by the engine, untouched, so that any
function of `edgeR`, `limma`, `variancePartition` or `DESeq2` can be
applied to it.

## Usage

``` r
fitObject(fit)

# S4 method for class 'RegionSetDE.fit'
fitObject(fit)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

## Value

The engine specific fit object.

## Author

Sebastian Gregoricchio
