# .linksFromGRanges

Turns a `GRanges` whose metadata columns name the genes into the long
link table the package stores. Convenient because this is the shape
`ChIPseeker::annotatePeak` and most contact callers return.

## Usage

``` r
.linksFromGRanges(linkRanges)
```

## Arguments

- linkRanges:

  `GRanges` with at least a gene identifier among its metadata columns.

## Value

A data.frame in the long link format.

## Author

Sebastian Gregoricchio
