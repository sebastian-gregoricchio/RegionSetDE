# .resolveLinkWeights

Decides what weights the regions of a gene, from what the links actually
carry. A weight measured or predicted for the link itself comes first,
the decay of the distance to the start site second, and nothing at all
last, which is a legitimate answer rather than a failure.

## Usage

``` r
.resolveLinkWeights(linkTable, decay = 5000)
```

## Arguments

- linkTable:

  Data.frame with the links of the regions kept.

- decay:

  Numeric value with the distance, in base pairs, at which the
  exponential decay falls to `1/e`.

## Value

A list with `weight`, a numeric vector parallel to `linkTable`, and
`source`, one among `"link.weight"`, `"distance"` and `"none"`.

## Author

Sebastian Gregoricchio
