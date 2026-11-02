# .asGeneSetList

Accepts the two shapes gene sets arrive in, a named list of gene
identifiers or a two-column table of the kind `clusterProfiler` calls a
`TERM2GENE`, and returns the list.

## Usage

``` r
.asGeneSetList(geneSets)
```

## Arguments

- geneSets:

  Named list of character vectors, or a data.frame with one column
  naming the set and one naming the gene.

## Value

A named list of character vectors.

## Author

Sebastian Gregoricchio
