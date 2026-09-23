# .bamWithIndex

Opens a BAM file with its index, handing the path over explicitly.
Rsamtools finds a BAI on its own and leaves a CSI alone, so a file
indexed with `samtools index -c` cannot be read by region unless it is
told where the index is. CSI is not an exotic case: BAI cannot address a
contig longer than 512 Mb at all, which rules it out for several plant
and amphibian assemblies.

## Usage

``` r
.bamWithIndex(bamFile)
```

## Arguments

- bamFile:

  String with the path of the BAM file.

## Value

A `BamFile` carrying the index that was found, or the path unchanged
when there is none.

## Author

Sebastian Gregoricchio
