# exportResults

Writes a result to disk: the table as a compressed TSV, the coordinates
as a BED that a genome browser will colour by direction, and every
parameter the analysis was run with as a flat file next to them.

## Usage

``` r
exportResults(
  results,
  path = ".",
  prefix = NULL,
  contrast = NULL,
  set = NULL,
  onlyChanging = FALSE,
  splitByDirection = FALSE,
  bedScore = "FDR",
  colourByStatus = TRUE,
  writeTiles = FALSE,
  provenance = TRUE,
  compress = TRUE,
  verbose = TRUE
)
```

## Arguments

- results:

  `RegionSetDE.results`, `RegionSetDE.setResults`, or either of the two
  list classes holding several contrasts.

- path:

  String with the directory the files are written to, created when it
  does not exist. Default: `"."`.

- prefix:

  String prepended to every file name. Default: `NULL`, the name of the
  contrast.

- contrast:

  String with the name of the contrast to write, or its position, when
  `results` holds several of them. Default: `NULL`, every one of them.

- set:

  Character vector with the names of the region sets to write. Default:
  `NULL`, all of them.

- onlyChanging:

  Logical value to indicate whether the files must hold only the regions
  labelled `"up"` or `"down"`. Default: `FALSE`.

- splitByDirection:

  Logical value to indicate whether a separate BED must be written for
  each direction, on top of the combined one. Default: `FALSE`.

- bedScore:

  String with the column mapped onto the BED score, one of `"FDR"`,
  `"log2FC"` and `"none"`. Default: `"FDR"`.

- colourByStatus:

  Logical value to indicate whether the BED must carry an item colour
  per direction, which makes it nine columns rather than six. Default:
  `TRUE`.

- writeTiles:

  Logical value to indicate whether the tile level table must be written
  as well, when the counts were tiled. Default: `FALSE`.

- provenance:

  Logical value to indicate whether the parameters of the analysis must
  be written alongside. Default: `TRUE`.

- compress:

  Logical value to indicate whether the tables must be gzipped. Default:
  `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

Invisibly, a character vector with the paths written.

## Details

The provenance file is the part worth keeping. It flattens everything
the object carries in its `parameters` slot into one parameter per line:
the counting level, the normalisation method and its factors, the filter
and its threshold, the engine, the design, the contrast, the multiple
testing correction, and on a fit with no replicates the dispersion and
where it came from. Six months later that file is the difference between
a result that can be reproduced and one that can only be repeated.

The BED is written with zero-based starts, as the format requires, so
its coordinates are one lower than the ones in the TSV. With
`colourByStatus` it carries nine columns and an item colour per
direction, which is what makes a browser show the increasing and
decreasing regions apart without a second file. The score maps
`-log10(FDR)` onto the range the format allows and saturates at the top
of it, so a score of 1000 means "at least that", not "exactly that".

A set level result has no coordinates of its own and only the table is
written.

## See also

[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
[`topRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/topRegions.md),
[`resultRanges`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultRanges.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

outputDirectory <- file.path(tempdir(), "regionSetDE-example")
dir.create(outputDirectory, showWarnings = FALSE)

exportResults(results, path = outputDirectory, prefix = "euratrans",
              verbose = FALSE)

list.files(outputDirectory)
#> [1] "euratrans_parameters.tsv" "euratrans_regions.bed"   
#> [3] "euratrans_regions.tsv.gz"

unlink(outputDirectory, recursive = TRUE)
```
