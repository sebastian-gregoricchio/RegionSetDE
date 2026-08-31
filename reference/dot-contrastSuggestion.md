# .contrastSuggestion

Builds the second half of the error message raised when a contrast
cannot be read, pointing at the reference level when that is what went
wrong.

## Usage

``` r
.contrastSuggestion(contrast, coefficientNames, colData = NULL)
```

## Arguments

- contrast:

  String with the contrast the user wrote.

- coefficientNames:

  Character vector with the columns of the design.

- colData:

  Sample metadata, or `NULL`.

## Value

A string, empty when nothing useful can be said.

## Author

Sebastian Gregoricchio
