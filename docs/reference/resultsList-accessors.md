# Accessors of RegionSetDE.resultsList

Extract one contrast, or the names of the contrasts, from a
`RegionSetDE.resultsList` object.

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
```

## Arguments

- x:

  `RegionSetDE.resultsList` object.

- i:

  String with the name of a contrast, or its position.

- j:

  Not used, present because the `[[` generic carries it.

- ...:

  Not used.

- name:

  String with the name of a contrast.

## Value

A `RegionSetDE.results` object, or a character vector for `names`.

## Author

Sebastian Gregoricchio
