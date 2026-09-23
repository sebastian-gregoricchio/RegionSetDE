# .correlationCellFunction

Builds the function writing the correlations in the cells of a heatmap,
in black or white depending on how dark the cell is.

## Usage

``` r
.correlationCellFunction(
  panelMatrix,
  digits = 2,
  fontSize = 9,
  valuesColour = NULL
)
```

## Arguments

- panelMatrix:

  Numeric matrix drawn by the heatmap.

- digits:

  Numeric value with the number of decimals written.

- fontSize:

  Numeric value with the font size of the names of the heatmap.

- valuesColour:

  String with a colour for every value, or `NULL`.

## Value

A function with the signature `ComplexHeatmap` expects for `cell_fun`.

## Author

Sebastian Gregoricchio
