# .plotSetScores

Draws a `RegionSetDE.setScores` object: one point per library per set,
the points of a library joined across the sets, and a bracket carrying
the paired difference and its adjusted p-value.

## Usage

``` r
.plotSetScores(
  setScores,
  set = NULL,
  groupBy = NULL,
  groupOrder = NULL,
  comparisons = NULL,
  annotate = TRUE,
  colours = NULL,
  baseColours = NULL,
  title = NULL,
  subtitle = NULL,
  legendPosition = "none",
  baseSize = 12
)
```

## Arguments

- setScores:

  `RegionSetDE.setScores` object.

- set:

  Character vector with the names of the region sets to draw. Default:
  `NULL`, all of them.

- groupBy:

  String with the name of a column of the sample metadata driving the
  colour. Default: `NULL`, the sample itself.

- groupOrder:

  Character vector with the levels of `groupBy` in the order they must
  appear. Default: `NULL`, alphabetical.

- comparisons:

  List of character vectors of length two, naming the pairs to annotate.
  Default: `NULL`, every pair that was compared.

- annotate:

  Logical value to indicate whether the brackets must be drawn. Default:
  `TRUE`.

- colours:

  Named character vector with one colour per sample. Default: `NULL`,
  shades built from the groups.

- baseColours:

  Character vector with one base colour per group. Default: `NULL`.

- title:

  String with the title of the plot, rendered as markdown. Default:
  `NULL`.

- subtitle:

  String with the subtitle of the plot, rendered as markdown. Default:
  `NULL`.

- legendPosition:

  String with the position of the legend. Default: `"none"`.

- baseSize:

  Numeric value with the base font size. Default: `12`.

## Value

A `ggplot` object.

## Details

The lines joining the points are not decoration. The comparison behind
the bracket is paired, one difference per library, and a figure drawing
the sets as independent clouds would be showing a test other than the
one whose p-value it carries. Three lines running parallel are what a
small difference with a decisive paired p-value looks like; three lines
crossing are why the interval on that difference is wide.

## Author

Sebastian Gregoricchio
