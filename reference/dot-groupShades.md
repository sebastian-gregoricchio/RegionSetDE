# .groupShades

Builds one colour per sample, a family per group and a shade per
replicate inside it.

## Usage

``` r
.groupShades(colTable, baseColours = NULL)
```

## Arguments

- colTable:

  Data.frame with the `sample` and `group` columns, already ordered.

- baseColours:

  Character vector with one base colour per group, or `NULL`.

## Value

A named character vector with one colour per sample.

## Author

Sebastian Gregoricchio
