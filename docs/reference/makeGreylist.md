# makeGreylist

Builds a greylist from input libraries: the stretches of genome where an
input carries far more fragments than the rest of its genome leads to
expect, such as copy number gains of the cell line, collapsed repeats or
regions that stick to any immunoprecipitation. Each input is judged
against its own coverage, and the regions flagged are pooled over the
inputs, ready for
[`applyGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyGreylist.md).

## Usage

``` r
makeGreylist(
  inputFiles,
  inputNames = NULL,
  binSize = 1024,
  quantile = 0.99,
  maxGap = 16384,
  minInputs = 1,
  excludeChromosomes = NULL,
  pairedEnd = "auto",
  fragmentLength = 150,
  maxFragmentLength = 1000,
  minMapq = 20,
  removeDuplicates = TRUE,
  nThreads = 1,
  verbose = TRUE
)
```

## Arguments

- inputFiles:

  Character vector with the paths of the input BAM files, or the
  data.frame returned by
  [`loadSampleSheet`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadSampleSheet.md),
  in which case every distinct file of its `input` column is used once
  and named after `input.id`.

- inputNames:

  Character vector with the names of the inputs. Default: `NULL`, the
  `input.id` of the sample sheet or the file names.

- binSize:

  Numeric value with the width of the windows, in base pairs.
  Consecutive windows overlap by half of it. Default: `1024`.

- quantile:

  Numeric value with the quantile of the fitted negative binomial above
  which a window is flagged. Default: `0.99`.

- maxGap:

  Numeric value, in base pairs: two flagged windows separated by a
  shorter gap are merged into one region. Default: `16384`.

- minInputs:

  Numeric value with the number of inputs that must flag a position for
  it to enter the greylist. Default: `1`, a position flagged by any
  input.

- excludeChromosomes:

  Character vector with the chromosomes left out, written in either
  naming style, as in
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).
  Default: `NULL`, none.

- pairedEnd:

  Logical value, one logical value per BAM file, or the string `"auto"`
  to read the layout from the files themselves. Default: `"auto"`.

- fragmentLength:

  Numeric value with the length to which single-end reads are extended.
  Default: `150`.

- maxFragmentLength:

  Numeric value with the maximum length accepted for a paired-end
  fragment. Default: `1000`.

- minMapq:

  Numeric value with the minimum mapping quality of a read. Default:
  `20`.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded. Default: `TRUE`.

- nThreads:

  Number of threads. Default: `1`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `GRanges` with one element per greylisted region, carrying `n.inputs`,
the number of inputs flagging it, and `inputs`, their names. Its
`metadata` holds `thresholds`, a data.frame with, for every input, the
fragments counted, the mean and size of the fitted negative binomial,
the threshold and how much of the genome was flagged, and `parameters`,
the arguments of the call.

## Details

The procedure follows GreyListChIP, which is also what DiffBind uses for
its greylists. The genome is cut into windows of `binSize` base pairs
overlapping by half, the fragments of each input are counted once in
each window holding their centre, and a negative binomial is fitted to
the counts of every input. Windows above the `quantile` of that
distribution are flagged, and flagged windows separated by less than
`maxGap` are merged into one region, the stretch between them included.

Two things differ from GreyListChIP. The distribution is fitted by
maximum likelihood on all the windows, rather than on 100 bootstrap
samples of 30,000 windows, so the threshold does not depend on the
random seed. And the fragments are read with the same filters as in
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md):
mapping quality, duplicates and proper pairs.

Empty windows enter the fit, as they do in GreyListChIP. A chromosome
without any read, such as chrY in a female sample or an unused contig,
adds a block of zeros that widens the fitted distribution and raises the
threshold, so it is worth leaving out through `excludeChromosomes`.

A greylist describes the input, not the chromatin. Marks like H3K9me3
sit largely on satellites and other repeats, which is precisely where
inputs pile up, so a greylist can take genuine signal away with the
artefacts. Check how much of each set
[`applyGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyGreylist.md)
removes, and consider `trimRegions = TRUE`, which cuts the greylisted
stretch out of a broad domain rather than dropping the whole domain.

## References

Brown G. GreyListChIP: Grey Lists – Mask Artefact Regions Based on ChIP
Inputs. Bioconductor package.

## See also

[`applyGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyGreylist.md),
[`loadSampleSheet`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadSampleSheet.md),
[`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
# The alignment shipped with Rsamtools stands in for an input library
inputFile <- system.file("extdata", "ex1.bam", package = "Rsamtools")

greylist <- makeGreylist(inputFile, binSize = 200, maxGap = 200)
#> Counting 1 input libraries over 32 windows of 200 bp...
#>   ex1: more than 331 fragments per window (mean 97.9), 0 windows flagged, 0 regions, 0 Mb.
#> Greylist: 0 regions, 0 Mb (0% of the genome), flagged by at least 1 of 1 inputs.
greylist
#> GRanges object with 0 ranges and 2 metadata columns:
#>    seqnames    ranges strand |  n.inputs      inputs
#>       <Rle> <IRanges>  <Rle> | <integer> <character>
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome

# The threshold and the fit behind it
S4Vectors::metadata(greylist)$thresholds
#>   input
#> 1   ex1
#>                                                                               file
#> 1 /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/Rsamtools/extdata/ex1.bam
#>   fragments   mean     size threshold windows flagged.windows regions
#> 1      1566 97.875 1.947997       331      32               0       0
#>   greylisted.bp
#> 1             0

if (FALSE) { # \dontrun{
# From a sample sheet, each distinct input once
sampleSheet <- loadSampleSheet("samples.csv")
greylist <- makeGreylist(sampleSheet, excludeChromosomes = c("chrM", "chrY"), nThreads = 4)

regions <- applyGreylist(regions, greylist = greylist)
} # }
```
