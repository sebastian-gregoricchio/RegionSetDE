# .registeredChromosomeNames

Converts chromosome names to the naming style of another set of names
through the tables of `GenomeInfoDb`, which cover the chromosomes of the
species it lists and their mitochondrion.

## Usage

``` r
.registeredChromosomeNames(chromosomeNames, targetSeqlevels)
```

## Arguments

- chromosomeNames:

  Character vector with the chromosome names to convert.

- targetSeqlevels:

  Character vector with chromosome names written in the style to convert
  to.

## Value

A character vector as long as `chromosomeNames`, with `NA` for the names
the tables do not hold, scaffolds and custom contigs for instance, or
when the style of `targetSeqlevels` is not recognised.

## Author

Sebastian Gregoricchio
