# RegionSetDE: differential chromatin signal over user-defined region sets

A peak caller defines the regions it then tests, and it redefines them
whenever the signal moves. A domain that spreads in the mutant comes
back as a wider peak, and the comparison quietly turns into a comparison
between two different sets of coordinates. `RegionSetDE` takes the
regions from the user instead. They are grouped into named sets of
arbitrary width and number, such as promoter classes, enhancer
catalogues, chromatin states or motif-derived sites, and they stay fixed
across conditions, so nothing about the regions depends on the outcome.

Most chromatin questions are asked about a class of elements rather than
about one locus. Whether H3K27ac drops at CpG island promoters is a
question about a set, and counting how many members of that set cross an
FDR threshold answers a different question: that count tracks sequencing
depth as much as it tracks the size of the effect. Here the set itself
is the unit of the test. One model fit serves both levels, so the
region-level and the set-level results can never disagree on the model
underneath them.

Set-level results come with a confidence interval on the condition
effect, computed from one score per library, which is where the
replication of the experiment lives, and next to it a second interval
for how much that effect varies from one region to the next. Thirty
thousand promoters are never counted as thirty thousand replicates.

## Entry points

- [`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md):

  reads BED, narrowPeak, broadPeak, `GRanges` or data.frames into named
  region sets.

- [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md),
  [`loadCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadCounts.md):

  count over the regions, or over fixed-width tiles of them, from BAM
  files, from bigWig coverage, or from a matrix that already exists.

- [`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md):

  scaling factors from background bins, from the regions themselves, or
  supplied from a spike-in or greenlist.

- [`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md):

  a single fit, with `edgeR`, `limma`-voom,
  [`variancePartition::dream`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md),
  `DESeq2` or `limma`-trend behind it.

- [`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md):

  tests one region at a time.

- [`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
  and
  [`testSetContrast`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md):

  test one whole set at a time, and the difference between two sets.

- [`makeSetUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md):

  builds the rows a set is compared against, matched on width and on
  baseline abundance.

- [`checkNullCalibration`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/checkNullCalibration.md):

  runs the contrast on rows taken to be null and shows whether their
  p-values come out uniform.

## Provenance

A substantial part of the code in this package was written with AI
assistance (Assisted-by: Claude, Anthropic). The methods, design and
validation are the author's, who maintains the package and is
responsible for its correctness.

## References

Wu D., Smyth G.K. (2012). Camera: a competitive gene set test accounting
for inter-gene correlation. *Nucleic Acids Research* 40(17):e133.

Lun A.T.L., Smyth G.K. (2016). csaw: a Bioconductor package for
differential binding analysis of ChIP-seq data using sliding windows.
*Nucleic Acids Research* 44(5):e45.

Robinson M.D., McCarthy D.J., Smyth G.K. (2010). edgeR: a Bioconductor
package for differential expression analysis of digital gene expression
data. *Bioinformatics* 26(1):139-140.

Law C.W., Chen Y., Shi W., Smyth G.K. (2014). voom: precision weights
unlock linear model analysis tools for RNA-seq read counts. *Genome
Biology* 15:R29.

## See also

Useful links:

- <https://github.com/sebastian-gregoricchio/RegionSetDE/>

- <https://sebastian-gregoricchio.github.io/RegionSetDE/>

- Report bugs at
  <https://github.com/sebastian-gregoricchio/RegionSetDE/issues>

## Author

Sebastian Gregoricchio
