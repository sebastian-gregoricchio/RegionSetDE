# .labelContrastTable

Puts the name of a contrast in front of its table, the way the tables of
several contrasts are stacked. The name is the one given to
`testRegions` and taken back by the `contrast` argument of every other
function, so that a stacked table can be filtered with the same string;
the description of the contrast follows it.

## Usage

``` r
.labelContrastTable(resultTable, contrastKey, singleResult)
```

## Arguments

- resultTable:

  Data.frame with the results of one contrast.

- contrastKey:

  String with the name of the contrast in the list.

- singleResult:

  The results object of that contrast.

## Value

The table with `contrast` and `contrast.description` as its first two
columns.

## Author

Sebastian Gregoricchio
