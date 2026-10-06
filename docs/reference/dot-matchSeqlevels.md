# .matchSeqlevels

Renames the chromosomes of a set of ranges so that they follow the names
of a signal file, or of another set of ranges. The chromosomes are
settled one by one: those already written as in the target are left
alone, and the others take the name the target gives them. Those the
target does not have under any name, scaffolds for the most part, follow
its style, so that the object does not come back half renamed. Only the
copy used for the counting is renamed, so the object returned to the
user keeps the style of the regions it was built from.

## Usage

``` r
.matchSeqlevels(x, targetSeqlevels, fileName = NULL, verbose = TRUE)
```

## Arguments

- x:

  `GRanges`, or any object accepting `seqlevels`, to be renamed.

- targetSeqlevels:

  Character vector with the chromosome names to align to, usually read
  from the header of a signal file.

- fileName:

  String with the file path, used in the messages. Default: `NULL`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

The input object with the renamed chromosomes.

## Author

Sebastian Gregoricchio
