# aggregateToGenes

Turns a per-region result into a per-gene table, through the links
stored by
[`annotateRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/annotateRegions.md).
One row per gene, with the regions that reached it counted in each
direction and the disagreement between them reported rather than
resolved.

This is the table a conventional over-representation analysis needs, and
the ranked list a gene-level GSEA needs, for the cases where the set
level route of
[`geneSetsToRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/geneSetsToRegionSets.md)
and
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
is not what is wanted. It is the weaker of the two routes and it is
honest about why.

## Usage

``` r
aggregateToGenes(
  results,
  method = "simes",
  decay = 5000,
  linkType = NULL,
  regionSets = NULL,
  adjustMethod = "BH",
  verbose = TRUE
)
```

## Arguments

- results:

  `RegionSetDE.results` object returned by
  [`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
  carrying links.

- method:

  Character vector with what to compute, among `"simes"`, which combines
  the p-values of the regions of a gene, and `"potential"`, which sums
  their fold changes weighted by how plausible their links are. Default:
  `"simes"`.

- decay:

  Numeric value with the distance, in base pairs, at which the weight of
  a link falls to `1/e`, used by `"potential"` when the links carry no
  weight of their own. Default: `5000`.

- linkType:

  Character vector with the kinds of link to follow, among `"promoter"`
  and `"distal"` as
  [`annotateRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/annotateRegions.md)
  labelled them with its `tssRegion`. This is how the promoter-proximal
  and the distal table are built separately. Default: `NULL`, every
  link.

- regionSets:

  Character vector with the names of the region sets the regions are
  taken from, which is the other way of restricting the table to one
  class of element, when the sets were already split at load time.
  Default: `NULL`, every set of the object.

- adjustMethod:

  String with the multiple testing correction applied across the genes.
  Default: `"BH"`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.geneResults` object. Its table, read with
[`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md),
holds one row per gene:

- `gene.id`, `gene.symbol`: the gene, and its symbol when the links
  carry one.

- `n.regions`: how many regions of the object reached it.

- `n.up`, `n.down`: how many of those are labelled `"up"` and `"down"`
  in the region table.

- `direction`: `"up"` or `"down"` when the regions that changed agree,
  `"mixed"` when they do not, `"unchanged"` when none of them changed.

- `log2FC`: the fold change of the region carrying the combined p-value,
  with `"simes"`.

- `p.value`, `FDR`: the Simes combination over the regions of the gene,
  corrected across the genes, with `"simes"`.

- `score`, `score.magnitude`: with `"potential"`, the weighted sum of
  the fold changes of the regions of the gene, `sum(w * log2FC)`, and
  the same sum on their absolute values, `sum(w * abs(log2FC))`.

## Details

The reason this function exists alongside
[`geneSetsToRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/geneSetsToRegionSets.md),
and the reason it is the second choice, is worth stating. Aggregating to
genes throws away the ranking of the regions and replaces several
measurements with one, which is exactly the step the set level route
avoids. Use it when an external tool needs a gene table, when a reviewer
asks for a conventional over-representation result, or when the gene
rather than the region is genuinely the unit of the question. For asking
which pathways responded, test the sets.

`n.up` and `n.down` are read from the `diff.status` column of the region
table rather than recomputed, so they honour the `FDR` and `log2FC` the
regions were labelled with, including after
[`updateThresholds`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/updateThresholds.md)
has changed them. The gene table therefore cannot disagree with the
region table about what counts as up.

`direction` is the answer to the question that makes this whole exercise
awkward. A gene reached by one region that opened and another that
closed is `"mixed"`, and nothing here decides which of the two wins,
because there is no defensible way to decide it: cancelling them loses a
real signal, and picking the stronger one is arbitrary. Treat the mixed
genes as their own category rather than dropping them. They are where a
rearranged regulatory landscape shows up, and a gene with three regions
up and three down is a different biological statement from a gene with
nothing happening, which a signed score cannot tell apart.

`"simes"` asks whether any region of the gene changed. It is relatively
relaxed, valid under the positive dependence that regions of one locus
have, and corrects for the number of regions a gene carries, which a
minimum p-value does not. The fold change reported beside it belongs to
the region that carried the p-value, not to an average, so a gene with
one strong region among many quiet ones is described by that region.

`"potential"` is a score and not a test: it has no error model and no
p-value of its own, which is why the Simes columns are worth keeping
alongside it. For a gene reached by regions `i` carrying weights `w[i]`,
`score` is `sum(w * log2FC)` and `score.magnitude` is
`sum(w * abs(log2FC))`. The first is signed, so a gene gaining at one
region and losing at another lands near zero; the second ignores sign
and does not, so the same gene scores high. Genes high in the second and
near zero in the first are the rearranged ones, the same signal
`direction` reports as a label.

`score` is the ranked gene-level statistic a GSEA takes, so
`fgsea::fgseaMultilevel(pathways = geneSets, stats = setNames(score, gene.id))`
is the whole of that route; ranking on `score.magnitude` asks which sets
were remodelled rather than which shifted, and
`sign(score) * -log10(p.value)` ranks on significance rather than on
effect. No wrapper for it is provided, since the call is one line and
wrapping it would add a dependency without adding a decision. What does
need stating in the methods is that the genes of the table are not
independent, so the permutation behind any such test is
anti-conservative.

The weights come from the links. A weight measured or predicted for the
link itself is used as it is, whatever it measures, which is what
`weightSource` in
[`annotateRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/annotateRegions.md)
records. Failing that, the distance to the start site is turned into
`exp(-d / decay)`, which is a guess about contact rather than a measure
of it. Failing both, every region of a gene counts the same. Which of
the three applied is in the `weighting` slot and is printed by `show`.

Two things this does not fix, both worth a line in the methods. Genes
are not independent here: a region linked to three genes contributes to
all three, so the rows of the table share data and any test run on them
afterwards is anti-conservative to an unquantified degree. And the
number of regions a gene carries depends on the length of its regulatory
locus, so a gene set rich in long-locus genes is easier to hit; a
count-based method built for that bias, `chipenrich` or `polyenrich`, is
the comparison to run when it matters.

A tiled result that was not combined back to regions is refused, since
its rows are tiles and a gene would then be weighted by how finely its
regions were cut. Test it with `combine = TRUE` first.

## See also

[`annotateRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/annotateRegions.md),
[`geneSetsToRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/geneSetsToRegionSets.md),
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md),
[`updateThresholds`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/updateThresholds.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
fit <- loadExampleData("fit", verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes

# Links of your own, four regions per gene
regionTable <- as.data.frame(SummarizedExperiment::rowRanges(fit@counts))[, 1:3]
geneNumber <- ceiling(nrow(regionTable) / 4)

annotatedFit <- annotateRegions(fit,
                                linkTable = data.frame(regionTable,
                                                       gene.id = rep(paste0("gene", seq_len(geneNumber)),
                                                                     each = 4)[seq_len(nrow(regionTable))]),
                                verbose = FALSE)

results <- testRegions(annotatedFit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

geneResults <- aggregateToGenes(results, method = c("simes", "potential"), verbose = FALSE)
geneResults
#> An object of class 'RegionSetDE.geneResults'
#>   contrast        : condition: SHR vs BN 
#>   engine          : edgeR 
#>   aggregation     : simes, potential 
#>   weighting       : none 
#>   genes           : 474 
#>   links           : 1895 
#>   direction       : up 3, down 9, mixed 0, unchanged 462 
#> 
#>  gene.id gene.symbol n.regions n.up n.down direction score  p.value      FDR
#>   gene54        <NA>         4    0      1      down -6.55 4.59e-08 2.18e-05
#>  gene166        <NA>         4    0      1      down -1.88 1.10e-05 2.60e-03
#>    gene4        <NA>         4    0      1      down -4.22 5.48e-05 6.58e-03
#>  gene309        <NA>         4    0      1      down -5.30 5.55e-05 6.58e-03
#>  gene295        <NA>         4    0      1      down -2.18 1.34e-04 1.27e-02
#>  gene162        <NA>         4    1      0        up  9.27 1.68e-04 1.33e-02

geneTable <- resultsTable(geneResults)
table(geneTable$direction)
#> 
#>      down unchanged        up 
#>         9       462         3 
```
