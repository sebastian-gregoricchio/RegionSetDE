# loadSampleSheet

Reads a table describing the libraries of an experiment, one row per
sample, and returns it with the file paths resolved and checked. Beside
the sample name and its signal file, a row can point to the peaks called
on that sample and to the input library it was sequenced against. Every
other column is kept as sample annotation, with no restriction on names
or number.

## Usage

``` r
loadSampleSheet(
  sampleSheet,
  columns = NULL,
  basePath = NULL,
  checkFiles = TRUE,
  verbose = TRUE
)
```

## Arguments

- sampleSheet:

  String with the path to a comma or tab separated file, or a
  data.frame.

- columns:

  Named character vector telling which column of the table holds each
  standard field, for instance `c(sample = "library", bam = "file")`.
  The standard fields are `sample`, `bam`, `bigwig`, `peaks`, `input`
  and `input.id`. A field not given here is looked for under its own
  name and then under the one used by DiffBind (`SampleID`, `bamReads`,
  `Peaks`, `bamControl`, `ControlID`), ignoring the case. Default:
  `NULL`.

- basePath:

  String with the directory the relative paths of the table refer to.
  Default: `NULL`, the directory holding the table when it is read from
  a file, the working directory otherwise.

- checkFiles:

  Logical value to indicate whether the existence of the files, and the
  index of the BAM files, must be checked. Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A data.frame with one row per sample. The standard fields found in the
table come first, under their standard names and with absolute paths,
followed by every other column as it was. Samples without an input or
without peaks carry `NA` there. When the table has inputs, `input.id`
names each distinct input library, taken from the table when it has such
a column and from the file names otherwise.

## Details

Only `sample` and one of `bam` and `bigwig` are required. The other
fields are optional, and so are their values: a sample can come without
peaks or without an input, and one input can serve several samples.
[`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md)
reads each distinct input file once, whatever the number of samples
pointing to it.

A DiffBind sample sheet is read as it is. `SampleID`, `bamReads`,
`Peaks`, `bamControl` and `ControlID` take the standard names, while
`Tissue`, `Factor`, `Condition`, `Treatment`, `Replicate` and any other
column stay as annotation, available to the design of
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
under their own names.

Relative paths are resolved against `basePath` and returned as absolute
paths, so the table stays valid when the working directory changes.
Addresses starting with a protocol, such as `https://`, are left
untouched and are not checked, since
[`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md)
can read remote bigWig files.

The table goes straight into the counting, as the source of the files,
of the sample names and of the annotation:
`countReads(regions, bamFiles = sampleSheet$bam, sampleNames = sampleSheet$sample, sampleMetadata = sampleSheet)`.

## See also

[`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md),
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
[`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
# Two samples sharing an input and a third one without any, paths relative to a project folder
sheetTable <- data.frame(sample = c("treated_1", "control_1", "control_2"),
                         bam = c("reads/treated_1.bam", "reads/control_1.bam", "reads/control_2.bam"),
                         peaks = c("peaks/treated_1.narrowPeak", "peaks/control_1.narrowPeak", NA),
                         input = c("reads/input_A.bam", "reads/input_A.bam", NA),
                         condition = c("treated", "control", "control"))

sampleSheet <- loadSampleSheet(sheetTable, basePath = "/data/project", checkFiles = FALSE)
#> Sample sheet with 3 samples. Peaks: 2 of 3. Inputs: 1 distinct, none for 1 of 3 samples.
sampleSheet
#>      sample                               bam
#> 1 treated_1 /data/project/reads/treated_1.bam
#> 2 control_1 /data/project/reads/control_1.bam
#> 3 control_2 /data/project/reads/control_2.bam
#>                                      peaks                           input
#> 1 /data/project/peaks/treated_1.narrowPeak /data/project/reads/input_A.bam
#> 2 /data/project/peaks/control_1.narrowPeak /data/project/reads/input_A.bam
#> 3                                     <NA>                            <NA>
#>   input.id condition
#> 1  input_A   treated
#> 2  input_A   control
#> 3     <NA>   control

# A DiffBind sample sheet is read as it is
diffbindTable <- data.frame(SampleID = c("MCF7_1", "MCF7_2"),
                            Condition = c("resistant", "responsive"),
                            Replicate = c(1, 1),
                            bamReads = c("reads/MCF7_1.bam", "reads/MCF7_2.bam"),
                            ControlID = c("MCF7_input", "MCF7_input"),
                            bamControl = c("reads/MCF7_input.bam", "reads/MCF7_input.bam"))

loadSampleSheet(diffbindTable, basePath = "/data/project", checkFiles = FALSE)
#> Columns read as standard fields: SampleID as sample, bamReads as bam, bamControl as input, ControlID as input.id.
#> Sample sheet with 2 samples. Peaks: none. Inputs: 1 distinct, none for 0 of 2 samples.
#>   sample                            bam                              input
#> 1 MCF7_1 /data/project/reads/MCF7_1.bam /data/project/reads/MCF7_input.bam
#> 2 MCF7_2 /data/project/reads/MCF7_2.bam /data/project/reads/MCF7_input.bam
#>     input.id  Condition Replicate
#> 1 MCF7_input  resistant         1
#> 2 MCF7_input responsive         1
```
