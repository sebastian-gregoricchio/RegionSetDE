# .matchChromosomeNames

Renames a plain vector of chromosome names into the naming style of a
signal file, the way `.matchSeqlevels` does for a set of ranges. What it
is for is `excludeChromosomes`: a name that matches nothing is not an
error, it simply excludes nothing, and since the chromosomes left out
decide the library sizes the silence would be paid for by the
normalisation.

## Usage

``` r
.matchChromosomeNames(
  chromosomeNames,
  targetSeqlevels,
  argumentName = "excludeChromosomes"
)
```

## Arguments

- chromosomeNames:

  Character vector with the chromosome names given by the user.

- targetSeqlevels:

  Character vector with the chromosome names of the files.

- argumentName:

  String with the name of the argument, used in the warning. Default:
  `"excludeChromosomes"`.

## Value

The names, converted where a conversion was needed.

## Author

Sebastian Gregoricchio
