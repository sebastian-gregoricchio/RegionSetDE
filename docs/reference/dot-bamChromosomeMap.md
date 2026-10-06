# .bamChromosomeMap

Brings the chromosomes of a group of BAM files under one set of names,
those of the first file or of a reference given by the caller. Files
aligned to the same assembly do not always name its chromosomes alike,
`chr1` in one header and `1` in the next, and the ranges being counted
can only be written one way. Every file is then read under its own names
and reported under the common ones.

## Usage

``` r
.bamChromosomeMap(
  bamFiles,
  referenceTargets = NULL,
  referenceLabel = NULL,
  verbose = FALSE
)
```

## Arguments

- bamFiles:

  Character vector with the paths of the BAM files.

- referenceTargets:

  Named vector with the chromosome lengths the names are taken from.
  Default: `NULL`, the header of the first file.

- referenceLabel:

  String naming the reference in the error messages. Default: `NULL`,
  the name of the first file.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `FALSE`.

## Value

A list with `lengths`, a named integer vector with the length of every
chromosome found in the reference or in a file, under the common names;
`names`, a list with one named character vector per file, giving the
name each of its chromosomes carries in the file (the names of the
vector are the common ones); `shared`, the common names of the
chromosomes every file has; and `renamed`, a logical vector telling
which files are read under names of their own.

## Details

Two files that give different lengths to the same chromosome are not on
the same assembly, and the counting stops there. A contig that some
files lack, under any name, is no reason to stop: scaffolds and decoys
differ between builds of the same assembly, and the files without it
simply hold no read there.

## Author

Sebastian Gregoricchio
