# .setAnnotation

Builds the label written above each region set, carrying the effect size
and the adjusted p-value of the set level test.

## Usage

``` r
.setAnnotation(setResults, setNames, yPosition)
```

## Arguments

- setResults:

  `RegionSetDE.setResults` object.

- setNames:

  Character vector with the sets present in the plot.

- yPosition:

  Numeric value with the height the labels are drawn at.

## Value

A data.frame with the `region.set`, `label` and `y.position` columns.

## Author

Sebastian Gregoricchio
