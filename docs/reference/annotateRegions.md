# annotateRegions

Links the regions of an object to genes and stores the links inside it.
The links are the weakest step of any functional analysis built on
chromatin data, so they are kept as data rather than resolved on the
fly: the full region-to-gene table travels with the object, one row per
link, and a region linked to three genes keeps all three.

Sources are accepted in three forms. Proximity to a transcription start
site, read either from a `TxDb` or an `EnsDb` through `txdb`, or from
ranges you already have through `genes`, which is what a GTF or a
species with no `TxDb` leaves you with. A table of links of your own
through `linkTable`, which is how the output of a contact-based method
arrives, Hi-C loops, ABC scores, universalEPI predictions or anything
else. Or a `GRanges` of regions carrying the gene identifiers in its
metadata columns, also through `linkTable`, which is the shape
[`ChIPseeker::annotatePeak`](https://rdrr.io/pkg/ChIPseeker/man/annotatePeak.html)
returns.

## Usage

``` r
annotateRegions(
  object,
  txdb = NULL,
  genes = NULL,
  orgdb = NULL,
  linkTable = NULL,
  maxDistance = NULL,
  flankDistance = NULL,
  tssRegion = c(-3000, 3000),
  keyType = "ENTREZID",
  weightSource = "supplied",
  verbose = TRUE
)
```

## Arguments

- object:

  `RegionSetDE.counts`, `RegionSetDE.fit`, `RegionSetDE.results` or
  `RegionSetDE.setResults` object. Not a `RegionSetDE` object: the links
  live beside the counted rows, so the regions have to be counted first.

- txdb:

  `TxDb` or `EnsDb` object with the gene models, used when the links are
  built by proximity. Default: `NULL`.

- genes:

  `GRanges` with the genes, or their start sites, carrying their
  identifier in a `gene.id` metadata column or in the names of the
  ranges, and optionally a `gene.symbol` column. The links are built by
  proximity to them, the way `txdb` does it, which is the route to take
  for an annotation read from a GTF or for a species with no `TxDb`.
  Default: `NULL`.

- orgdb:

  `OrgDb` object used to add the gene symbols. Default: `NULL`, no
  symbols.

- linkTable:

  Data.frame or `GRanges` with links of your own. Default: `NULL`, the
  links are built from `txdb`.

- maxDistance:

  Numeric value with the distance beyond which no gene is assigned, in
  base pairs, when the nearest gene is taken. Regions further than this
  are recorded with a missing `gene.id` rather than dropped. Default:
  `NULL`, every region takes its nearest site whatever the distance.

- flankDistance:

  Numeric value: every gene whose start site falls within it is linked
  to the region, rather than the nearest one alone. This is what makes
  the proximity route many-to-many, and it is how a region over a
  bidirectional promoter keeps both of its genes. Regions reaching
  nothing are recorded unassigned. Cannot be used with `maxDistance`,
  which caps the nearest-gene route instead. Default: `NULL`, the
  nearest site.

- tssRegion:

  Numeric vector of length two with the window around a start site,
  signed and read on the strand of the gene, inside which a link is
  called a promoter one. It fills the `link.type` column and decides
  nothing else. Default: `c(-3000, 3000)`, the window `ChIPseeker` uses.

- keyType:

  String with the key type the gene identifiers follow in `orgdb`.
  Default: `"ENTREZID"`.

- weightSource:

  String naming where the weights of `linkTable` come from, stored next
  to them. Default: `"supplied"`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

The input object with the links stored in it. The full table is in
`metadata(x)$links` and is read with
[`regionLinks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/regionLinks.md),
while `gene.id`, `gene.symbol`, `distance.to.tss` and `n.links` are
written next to the regions as one representative link per region, for
reading and for labelling plots.

## Details

A region-to-gene link is a prediction, and the quality of every
gene-level result rests on it rather than on the statistics applied
afterwards. Nearest transcription start site is wrong for a substantial
share of distal elements, which skip over genes and reach promoters
megabases away. That is the reason this function stores a table instead
of a column: the proximity links can be replaced by measured or
predicted contacts without recomputing anything else, and the two can be
compared on the same object.

The links can be many-to-many, and whether they are is a choice. A
supplied `linkTable` carries whatever links it was built with, so a
region overlapping a bidirectional promoter keeps both genes and a
region inside an intron keeps the host gene alongside whatever else its
contacts reach. The proximity route takes the single nearest start site
by default, one link per region, and `flankDistance` is what turns it
many-to-many: every gene whose start site falls within that distance is
linked, each with its own signed distance. Either way nothing decides
which link is right, and the number of links is not the number of
regions nor the number of genes.

`tssRegion` fills `link.type` with `"promoter"` or `"distal"`, the
window being read on the strand of the gene so that a region upstream of
a minus-strand gene sits at a higher coordinate than its start site. The
column decides nothing by itself; it is there so that
[`geneSetsToRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/geneSetsToRegionSets.md)
and
[`aggregateToGenes`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/aggregateToGenes.md)
can run the promoter-proximal and the distal analysis separately, which
is usually worth doing, since much of the apparent disagreement within a
gene is promoters staying put while enhancers move. What this function
does not do is classify a region into a genomic feature, promoter
against intron against exon, with a priority order; that is what
[`ChIPseeker::annotatePeak`](https://rdrr.io/pkg/ChIPseeker/man/annotatePeak.html)
is for, and its output goes into `linkTable` directly.

Weights are optional and they carry no assumed scale. An ABC score, a
Hi-C loop score, a predicted contact probability and a distance decay
are all accepted as they are, with `weightSource` recording which of
them it was, and nothing is normalised. Supplying no weight at all is
the normal case and is fully supported: the links are then unweighted,
and every function reading them treats the links of a region as equally
plausible. What a weight cannot do is reach the set-level tests, since
neither `camera` nor `fry` accepts a per-region weight. Where it does
reach is
[`geneSetsToRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/geneSetsToRegionSets.md),
which uses it as a confidence filter through its `minWeight` argument.

Two kinds of mismatch between the table and the object are handled
differently on purpose. Links on regions the object does not hold are
dropped, because a contact map is routinely computed over a wider peak
set than the one being analysed. Regions the table says nothing about
are kept, with a missing gene, because dropping them would shrink the
universe of every enrichment run afterwards without anything in the
output saying so. The counts of both are reported and recorded in the
parameters.

`maxDistance` makes a third case explicit. A region whose nearest start
site is 800 kb away has a nearest gene in the arithmetic sense only, and
assigning it is worse than leaving it alone. Those regions are written
into the table with a missing `gene.id`, which keeps them in the
universe, keeps them out of every gene set, and leaves the decision
visible in the object rather than buried in the choice of which rows
were loaded.

Gene identifiers are matched exactly, as strings, and must follow the
same namespace as the gene sets they will be tested against. Symbols are
stored when an `OrgDb` is supplied and are never used for matching: they
are not stable between annotation releases and not unique within one.

## See also

[`regionLinks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/regionLinks.md),
[`geneSetsToRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/geneSetsToRegionSets.md),
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)

# Links of your own, the shape a contact-based method returns. Weights are optional.
regionTable <- as.data.frame(SummarizedExperiment::rowRanges(counts))[1:6, 1:3]
linkTable <- data.frame(regionTable,
                        gene.id = c("11287", "11287", "66968", "73824", "73824", "20683"),
                        link.weight = c(0.42, 0.08, 0.31, 0.57, 0.12, 0.66))

counts <- annotateRegions(counts, linkTable = linkTable,
                          weightSource = "ABC score", verbose = FALSE)

regionLinks(counts)
#>                                        region gene.id link.weight weight.source
#> promoterNonCpG|region_00002   chr12:4041-5040   11287        0.42     ABC score
#> promoterNonCpG|region_00003   chr12:6212-7211   11287        0.08     ABC score
#> promoterNonCpG|region_00005 chr12:10685-11684   66968        0.31     ABC score
#> promoterNonCpG|region_00006 chr12:13540-14539   73824        0.57     ABC score
#> promoterNonCpG|region_00008 chr12:18501-19500   73824        0.12     ABC score
#> promoterNonCpG|region_00009 chr12:20855-21854   20683        0.66     ABC score
head(SummarizedExperiment::rowData(counts)$gene.id)
#> [1] "11287" "11287" "66968" "73824" "73824" "20683"
```
