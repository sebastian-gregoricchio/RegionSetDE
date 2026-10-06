# .translateChromosomeNames

Finds, for every chromosome name, the name the same chromosome goes
under in another set of names. A name already there is kept. For the
others the conversions of `GenomeInfoDb` are tried first, then the `chr`
prefix is added or dropped, and the mitochondrion is looked for under
the four names it is given (`chrM`, `MT`, `chrMT`, `M`). Each name is
settled on its own, so a set mixing two styles is translated as well as
a set written in one.

## Usage

``` r
.translateChromosomeNames(chromosomeNames, targetSeqlevels)
```

## Arguments

- chromosomeNames:

  Character vector with the chromosome names to translate.

- targetSeqlevels:

  Character vector with the chromosome names to translate into.

## Value

A character vector as long as `chromosomeNames`, with the matching name
of `targetSeqlevels` and `NA` where there is none.

## Author

Sebastian Gregoricchio
