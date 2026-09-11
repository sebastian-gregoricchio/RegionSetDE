# .resolveContrast

Turns the `contrast` argument of
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
into a numeric vector over the columns of the design.

## Usage

``` r
.resolveContrast(contrast, design, colData = NULL)
```

## Arguments

- contrast:

  String with a coefficient name or an expression over the coefficients,
  a character vector of length three naming a column and two of its
  levels, or a numeric vector.

- design:

  Design matrix.

- colData:

  `DataFrame` or data.frame with the sample metadata, needed by the
  three-element form. Default: `NULL`.

## Value

A list with the `vector` of coefficients, a `label` describing the
contrast, and the `column` of the `colData` and the two `groups` it
separates when the contrast turns out to be a difference between two
levels of one variable.

## Author

Sebastian Gregoricchio
