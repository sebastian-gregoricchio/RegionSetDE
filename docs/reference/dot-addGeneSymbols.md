# .addGeneSymbols

Adds the gene symbols to a link table, read from an `OrgDb`. Cosmetic:
the symbols are for reading and for labelling, nothing is ever matched
on them, since they are neither stable between releases nor unique.

## Usage

``` r
.addGeneSymbols(linkTable, orgdb, keyType = "ENTREZID", verbose = TRUE)
```

## Arguments

- linkTable:

  Data.frame with the links.

- orgdb:

  `OrgDb` object.

- keyType:

  String with the key type of `gene.id` in the `OrgDb`. Default:
  `"ENTREZID"`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

The link table with a `gene.symbol` column.

## Author

Sebastian Gregoricchio
