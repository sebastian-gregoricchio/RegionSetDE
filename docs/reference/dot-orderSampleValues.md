# .orderSampleValues

Checks a vector of user supplied values against the samples of a counts
object and puts it in the order of the columns.

## Usage

``` r
.orderSampleValues(values, sampleNames, parameterName)
```

## Arguments

- values:

  Numeric vector, named or not.

- sampleNames:

  Character vector with the sample names, in the order of the columns.

- parameterName:

  String with the name of the parameter, used in the error messages.

## Value

The numeric vector, ordered as the columns of the object.

## Author

Sebastian Gregoricchio
