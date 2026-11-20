# Accessors of RegionSetDE.resultsList

Extract one contrast, all of them, or the names of the contrasts, from a
`RegionSetDE.resultsList` or a `RegionSetDE.setResultsList` object.

## Usage

``` r
# S4 method for class 'RegionSetDE.setResultsList,ANY,ANY'
x[[i, j, ...]]

# S4 method for class 'RegionSetDE.setResultsList'
x$name

# S4 method for class 'RegionSetDE.setResultsList'
names(x)

# S4 method for class 'RegionSetDE.setResultsList'
length(x)

# S4 method for class 'RegionSetDE.resultsList,ANY,ANY'
x[[i, j, ...]]

# S4 method for class 'RegionSetDE.resultsList'
x$name

# S4 method for class 'RegionSetDE.resultsList'
names(x)

# S4 method for class 'RegionSetDE.resultsList'
length(x)

# S4 method for class 'RegionSetDE.resultsList'
as.list(x, ...)

# S4 method for class 'RegionSetDE.setResultsList'
as.list(x, ...)
```

## Arguments

- x:

  `RegionSetDE.resultsList` or `RegionSetDE.setResultsList` object.

- i:

  String with the name of a contrast, or its position.

- j:

  Not used, present because the `[[` generic carries it.

- ...:

  Not used.

- name:

  String with the name of a contrast.

## Value

A `RegionSetDE.results` object, or a `RegionSetDE.setResults` object for
the list of set level results. `names` returns a character vector,
`length` the number of contrasts, and `as.list` a named list with every
result, one per contrast.

## Author

Sebastian Gregoricchio
