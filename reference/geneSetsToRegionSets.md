# geneSetsToRegionSets

Turns gene sets into region sets, by mapping every gene of a set through
the links stored by
[`annotateRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/annotateRegions.md),
so that a pathway becomes the regions assigned to its genes and can be
handed straight to
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md).

This is what makes pathway-level testing of differential regions go
through `camera` and `fry` rather than through a ranked gene list. The
directionality problem that ordinarily has to be solved before a gene
list can be written, one gene carrying several regions moving in
different directions, is not solved here: it is left to the test, which
reports a self-contained and a mixed p-value and so answers it rather
than deciding it in advance.

## Usage

``` r
geneSetsToRegionSets(
  object,
  geneSets,
  minWeight = NULL,
  linkType = NULL,
  minRegions = 2,
  maxRegions = NULL,
  regionSets = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  `RegionSetDE.counts`, `RegionSetDE.fit` or `RegionSetDE.results`
  object carrying the links.

- geneSets:

  Named list of gene identifiers, one element per set, or a two-column
  table naming the set and the gene.

- minWeight:

  Numeric value with the smallest `link.weight` a link must carry to be
  used. Default: `NULL`, every link is used whatever its weight.

- linkType:

  Character vector with the kinds of link to follow, among `"promoter"`
  and `"distal"` as
  [`annotateRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/annotateRegions.md)
  labelled them with its `tssRegion`. This is how the promoter-proximal
  and the distal analysis are run separately. Default: `NULL`, every
  link.

- minRegions:

  Numeric value with the smallest number of regions a set must reach to
  be returned. Note that this counts regions and not genes. Default:
  `2`, the fewest a set level test can work with.

- maxRegions:

  Numeric value with the largest number of regions a set may reach.
  Default: `NULL`, no upper bound.

- regionSets:

  Character vector with the names of the region sets the regions are
  taken from, which is how the analysis is restricted to one class of
  element, promoters for instance. Default: `NULL`, every set of the
  object.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A named list of character vectors, one element per gene set, holding the
keys of the regions reached by that set. The keys are the set name and
the region identifier pasted with a vertical bar, which is the form
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
accepts through its `regionSets` argument. A `geneSetSummary` attribute
carries one row per set with the number of genes it holds, the number of
them the links reach, and the number of regions.

## Details

Gene sets overlap and region sets of this package do not. A promoter
belonging to fifty pathways is one row of the object that fifty sets
reach, which is why this returns keys rather than a `GRangesList`: a
list of ranges would have to repeat the row fifty times, and the
counting would then be wrong in fifty ways.
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
resolves the keys against its own rows, after any collapsing of the
tiles, so the sets stay correct whatever the counting level.

The comparison universe follows from that. With sets supplied this way,
each one is compared against every other row of the object rather than
against the other sets, which is the right comparison for a pathway: the
question is whether the regions of this pathway moved more than the
regions of the rest of the genome that was tested, not whether they
moved more than the regions of the other pathways. It also means the
universe is the annotated and tested regions and nothing wider, which is
the background any enrichment of this kind should have and the one most
often got wrong.

What this does not fix is the bias that comes before it. A gene with a
long regulatory locus carries more regions than a gene with a short one,
so a set rich in such genes reaches more rows, and `camera` corrects for
the correlation between those rows without correcting for how many of
them there are. Where that matters, compare the result against a
count-based method built for it, `chipenrich` or `polyenrich` for
instance, run on the regions rather than on the genes.

`minWeight` is the one place a link weight is used, and it is used as a
filter rather than as a weight: neither `camera` nor `fry` accepts a
per-region weight, so a contact probability can decide whether a link is
trusted enough to be followed and cannot make one link count more than
another. Asking for it on links that carry no weight stops, rather than
silently keeping everything, since a filter that quietly does nothing is
worse than one that is not there.

Splitting promoter-proximal from distal links is usually worth doing,
and there are two ways to reach it. `linkType` filters the links
themselves, on the labels
[`annotateRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/annotateRegions.md)
wrote from its `tssRegion`, which is the route when the regions are one
undifferentiated set of peaks. `regionSets` selects whole region sets
instead, which is the route when they were already split at load time,
by `loadConsensusPeaks(regionSets = , regionMode = "split")` for
instance. Either way the point is the same: promoter-proximal and distal
regions of one gene behave differently often enough that mixing them is
where much of the apparent disagreement within a gene comes from.

## See also

[`annotateRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/annotateRegions.md),
[`regionLinks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/regionLinks.md),
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)

regionTable <- as.data.frame(SummarizedExperiment::rowRanges(counts))[1:8, 1:3]
linkTable <- data.frame(regionTable,
                        gene.id = c("geneA", "geneA", "geneB", "geneB",
                                    "geneC", "geneC", "geneD", "geneD"))

counts <- annotateRegions(counts, linkTable = linkTable, verbose = FALSE)

geneSets <- list(pathwayOne = c("geneA", "geneB"),
                 pathwayTwo = c("geneC", "geneD", "geneE"))

regionSets <- geneSetsToRegionSets(counts, geneSets = geneSets, verbose = FALSE)
lengths(regionSets)
#> pathwayOne pathwayTwo 
#>          4          4 
attr(regionSets, "geneSetSummary")
#>     set.name n.genes n.genes.linked n.regions
#> 1 pathwayOne       2              2         4
#> 2 pathwayTwo       3              2         4
```
