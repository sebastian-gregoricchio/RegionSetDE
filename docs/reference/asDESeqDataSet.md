# asDESeqDataSet

Turns a `RegionSetDE.counts` object into a `DESeqDataSet`, with the
normalisation stored in the object carried across as normalisation
factors.

## Usage

``` r
asDESeqDataSet(
  counts,
  design = NULL,
  assay = "counts",
  useOffsets = TRUE,
  verbose = TRUE
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object, or a `RegionSetDE.fit`.

- design:

  Formula evaluated on the `colData`, the same formula written as a
  string, or a design matrix. Default: `NULL`, the design of the fit
  when one is given and `~ 1` otherwise.

- assay:

  String with the name of the assay holding the counts. Default:
  `"counts"`.

- useOffsets:

  Logical value to indicate whether the normalisation stored in the
  object must be carried across. Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `DESeqDataSet`, unfitted.

## Details

The normalisation goes into `normalizationFactors` rather than
`sizeFactors`, since the object may hold one value per region and per
sample rather than one per sample. `DESeq2` asks those factors to have a
geometric mean of one on each row, which is what the row centring here
is for; without it the fitted values come out on a scale that has
nothing to do with the counts. The object is returned unfitted, so
[`DESeq2::DESeq`](https://rdrr.io/pkg/DESeq2/man/DESeq.html) still has
to be run on it. Nothing stops a design being passed that differs from
the one the package used, which is the point of the conversion, but a
result obtained that way is a different analysis and not a check of this
one. There is no `as(counts, "DESeqDataSet")` to go with it, because a
coercion has to be registered when the package is built and the class it
points at only exists when `DESeq2` is installed.

## See also

[`asDGEList`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/asDGEList.md),
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes

if (requireNamespace("DESeq2", quietly = TRUE)) {
  deseqData <- asDESeqDataSet(counts, design = ~ condition, verbose = FALSE)
  deseqData
}
#> class: DESeqDataSet 
#> dim: 3224 4 
#> metadata(1): version
#> assays(2): counts normalizationFactors
#> rownames(3224): promoterNonCpG|region_00002 promoterNonCpG|region_00003
#>   ... promoterCpG|region_03797 promoterCpG|region_03798
#> rowData names(3): region.set region.id tile.id
#> colnames(4): lv-H3K4me3-BN-female-bio1-tech1
#>   lv-H3K4me3-BN-male-bio2-tech1 lv-H3K4me3-SHR-male-bio2-tech1
#>   lv-H3K4me3-SHR-male-bio3-tech1
#> colData names(9): sample bam.file ... norm.factor scaling.factor
```
