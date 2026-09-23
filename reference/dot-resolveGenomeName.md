# .resolveGenomeName

Turns the name of an assembly into the one the files are indexed under,
so that GRCh38 and hg38 reach the same list.

## Usage

``` r
.resolveGenomeName(genome)
```

## Arguments

- genome:

  String with the assembly given by the user.

## Value

A string with the name used in the index, the input itself when it is
not a known alias.

## Author

Sebastian Gregoricchio
