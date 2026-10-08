# updateObject methods for the RegionSetDE classes

Bring an object saved by an earlier version of RegionSetDE to the
current definition of its class, as
[`BiocGenerics::updateObject`](https://rdrr.io/pkg/BiocGenerics/man/updateObject.html)
does for the Bioconductor classes. An object saved before the `greylist`
slot existed gets it, filled with the greylist of the consensus when the
regions were built by
[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md),
and empty otherwise. A fit or a result updates the counts it carries as
well, and a list of results every result in it.

## Usage

``` r
# S4 method for class 'RegionSetDE'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'RegionSetDE.counts'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'RegionSetDE.fit'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'RegionSetDE.results'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'RegionSetDE.setResults'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'RegionSetDE.setScores'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'RegionSetDE.resultsList'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'RegionSetDE.setResultsList'
updateObject(object, ..., verbose = FALSE)
```

## Arguments

- object:

  An object of one of the RegionSetDE classes.

- ...:

  Passed on to the method of `SummarizedExperiment` for a counts object,
  ignored otherwise.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `FALSE`.

## Value

The object, updated to the current definition of its class. An object
already up to date comes back as it was.

## Details

A greylist applied with
[`applyGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyGreylist.md)
before the slot existed left only its size in the object, so the update
cannot recover it: apply it again, or pass it to the `greylist` argument
of
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
counts <- BiocGenerics::updateObject(counts, verbose = TRUE)
#> updateObject(object="ANY") default for object of class 'matrix'
#> [updateObject] Validating the updated object ... 
#> OK
#> [updateObject] Validating the updated object ... 
#> OK
#> [updateObject] DFrame object is current.
#> [updateObject] Nothing to update.
#> updateObject(object="ANY") default for object of class 'NULL'
#> [updateObject] Validating the updated object ... 
#> OK
#> [updateObject] Validating the updated object ... 
#> OK
#> [updateObject] DFrame object is current.
#> [updateObject] Nothing to update.
#> updateObject(object="ANY") default for object of class 'NULL'
#> [updateObject] Validating the updated object ... 
#> OK
#> [updateObject] Validating the updated object ... 
#> OK
#> [updateObject] Internal representation of GRanges object is current.
#> [updateObject] Nothing to update.
#> [updateObject] Internal representation of IRanges object is current.
#> [updateObject] Nothing to update.
#> updateObject(object="ANY") default for object of class 'NULL'
#> [updateObject] Validating the updated object ... 
#> OK
#> [updateObject] Validating the updated object ... 
#> OK
#> [updateObject] DFrame object is current.
#> [updateObject] Nothing to update.
#> updateObject(object="ANY") default for object of class 'NULL'
#> [updateObject] Validating the updated object ... 
#> OK
#> [updateObject] Validating the updated object ... 
#> OK
#> [updateObject] Validating the updated object ... 
#> OK
#> [updateObject] Validating the updated object ... 
#> OK
```
