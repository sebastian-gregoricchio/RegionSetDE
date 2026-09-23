# .bamIndexPath

Looks for the index of a BAM file under the four names one can go under:
`file.bam.bai`, `file.bai`, `file.bam.csi` and `file.csi`.

## Usage

``` r
.bamIndexPath(bamFile)
```

## Arguments

- bamFile:

  String with the path of the BAM file.

## Value

String with the path of the first index found, `NA` when the file has
none.

## Author

Sebastian Gregoricchio
