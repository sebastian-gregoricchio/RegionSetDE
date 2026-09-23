# .styleSeqlevels

Renames the chromosomes of a `GRanges` to a naming style, editing the
prefix by hand when the conversion of `GenomeInfoDb` does not cover the
contigs at hand.

## Usage

``` r
.styleSeqlevels(x, seqlevelsStyle = "UCSC")
```

## Arguments

- x:

  `GRanges` object.

- seqlevelsStyle:

  String with the style, one among `"UCSC"`, `"Ensembl"` and `"NCBI"`.

## Value

The object with its chromosomes renamed.

## Author

Sebastian Gregoricchio
