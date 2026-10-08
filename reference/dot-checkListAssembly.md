# .checkListAssembly

Refuses a list of regions built for another assembly than the one
declared for the peaks, which would overlap them on the chromosome names
alone.

## Usage

``` r
.checkListAssembly(listInput, genomeAssembly, listLabel, regionLabel = "peaks")
```

## Arguments

- listInput:

  The list as given: a `GRanges`, a path, a data.frame, or a list of
  them.

- genomeAssembly:

  String with the assembly of the peaks, or `NULL`.

- listLabel:

  String with the name of the parameter the list came from.

- regionLabel:

  String naming what the list is applied to in the message. Default:
  `"peaks"`.

## Value

Nothing, it stops when the two assemblies differ.

## Author

Sebastian Gregoricchio
