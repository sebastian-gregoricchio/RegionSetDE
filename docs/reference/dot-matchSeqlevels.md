# .matchSeqlevels

Renames the chromosomes of a set of ranges so that they follow the
naming style of a signal file, leaving them untouched when the two
already agree. Only the copy used for the counting is renamed, so the
object returned to the user keeps the style of the regions it was built
from.

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
