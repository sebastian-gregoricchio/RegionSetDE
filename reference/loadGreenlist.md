# loadGreenlist

Returns a CUT&RUN or CUT&Tag greenlist shipped with the package, the
regions whose background is consistent enough between experiments to
normalise on. It goes to
[`countGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countGreenlist.md),
which counts the libraries over it before
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)
turns those counts into scaling factors.

## Usage

``` r
loadGreenlist(
  genome,
  assay,
  source = NULL,
  seqlevelsStyle = "UCSC",
  verbose = TRUE
)
```

## Arguments

- genome:

  String with the genome assembly, `"hg38"` or `"mm39"` for the
  published lists. The usual aliases are understood, `"GRCh38"` for
  instance.

- assay:

  String with the assay the list was built for, `"cutrun"` or
  `"cuttag"`.

- source:

  String with the laboratory or project the list comes from, as
  [`availableRegionLists`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/availableRegionLists.md)
  prints it, needed only when a genome carries more than one list for
  the same assay. Default: `NULL`.

- seqlevelsStyle:

  String with the chromosome naming style of the output, one among
  `"UCSC"` (chr1), `"Ensembl"` (1) and `"NCBI"`, or `NULL` to keep the
  names of the file. Default: `"UCSC"`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `GRanges` with the regions of the greenlist and their name, carrying
the assembly in its `genome` and the source, the version and the
reference in its metadata.

## Details

The greenlist was built by measuring, over hundreds of public libraries,
which 1 kb bins carry background of a consistent magnitude, keeping the
most consistent of them and discarding anything within 5 kb of a gene so
that genuine signal is not called noise. The reads landing there follow
the amount of material sequenced rather than the factor being mapped,
which is what makes them usable as an internal reference when no
spike-in was added.

The list is specific to the protocol. CUT&RUN and CUT&Tag have their
own, and using one for the other means normalising on regions whose
background was never shown to be consistent in that assay. There is no
greenlist for mm10, only for mm39, so an mm10 analysis needs the mm39
list lifted over, and the lifted file should be kept beside the chain it
came from rather than passed off as the published one.

## See also

[`countGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countGreenlist.md),
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md),
[`loadBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadBlacklist.md),
[`availableRegionLists`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/availableRegionLists.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
greenlist <- loadGreenlist("hg38", assay = "cutrun")
#> The deMello greenlist v1 for hg38: 869 regions covering 1.7 Mb.
length(greenlist)
#> [1] 869
sum(BiocGenerics::width(greenlist))
#> [1] 1725126
```
