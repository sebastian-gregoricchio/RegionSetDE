# .checkGreenlistLibraries

Warns about the libraries whose greenlist counts are too thin to give a
steady factor: those below the depth the lists were built from, and
those covering far fewer regions of the list than the others.

## Usage

``` r
.checkGreenlistLibraries(
  greenlistExperiment,
  librarySizes,
  greenlistInfo = NULL
)
```

## Arguments

- greenlistExperiment:

  `SummarizedExperiment` with the greenlist counts, as built by
  `countGreenlist`.

- librarySizes:

  Numeric vector with the library size of every sample, possibly `NA`.

- greenlistInfo:

  List with the metadata of the greenlist, as set by
  [`loadGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadGreenlist.md).

## Value

Nothing, it raises warnings at most.

## Author

Sebastian Gregoricchio
