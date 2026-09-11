# .asFormula

Parses a design written as a string, leaving anything else untouched.

## Usage

``` r
.asFormula(x, parameterName = "design")
```

## Arguments

- x:

  Object passed by the user.

- parameterName:

  String with the name of the parameter, used in the error message.

## Value

A formula when the input was a string, the input itself otherwise.

## Author

Sebastian Gregoricchio
