# .sampleIndex

Turns a selection of samples, given by name, by position or as a logical
vector, into the positions of the samples kept, in the order of the
object.

## Usage

``` r
.sampleIndex(sampleNames, samples = NULL)
```

## Arguments

- sampleNames:

  Character vector with the names of all the samples, in the order of
  the object.

- samples:

  Character, numeric or logical vector with the samples kept, or `NULL`
  for all of them.

## Value

An integer vector with the positions of the samples kept.

## Author

Sebastian Gregoricchio
