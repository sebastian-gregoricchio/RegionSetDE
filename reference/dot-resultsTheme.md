# .resultsTheme

Extends the theme shared by the package with the pieces the result plots
need: a centred subtitle rendered as markdown, axis labels in black at
full size, and optionally angled labels on the x axis. The weight of the
subtitle is written out rather than left to the inheritance, which
otherwise picks up the bold of the title.

## Usage

``` r
.resultsTheme(legendPosition = "right", baseSize = 12, rotateX = FALSE)
```

## Arguments

- legendPosition:

  String with the position of the legend. Default: `"right"`.

- baseSize:

  Numeric value with the base font size. Default: `12`.

- rotateX:

  Logical value to indicate whether the labels of the x axis must be
  angled. Default: `FALSE`.

## Value

A `ggplot2` theme.

## Author

Sebastian Gregoricchio
