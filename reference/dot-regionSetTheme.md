# .regionSetTheme

Theme shared by the plotting functions of the package. It follows the
look of a publication panel: white background, no grid, black axis
lines, no box around the facet labels. Titles and axis titles are
rendered as markdown, so that subscripts such as `log<sub>2</sub>` come
out formatted.

## Usage

``` r
.regionSetTheme(legendPosition = "right", baseSize = 12)
```

## Arguments

- legendPosition:

  String with the position of the legend, one among `"right"`, `"left"`,
  `"top"`, `"bottom"` or `"none"`. Default: `"right"`.

- baseSize:

  Numeric value with the base font size. Default: `12`.

## Value

A `ggplot2` theme.

## Author

Sebastian Gregoricchio
