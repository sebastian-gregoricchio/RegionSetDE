# .validateLinkTable

Checks a region-to-gene link table and brings it onto the chromosome
names of the regions. The two directions of mismatch are not treated
alike: links landing on regions the object does not hold are dropped,
since a contact map is routinely computed on a wider set of peaks than
the one being tested, while regions the table says nothing about are
kept with a missing gene, since dropping them would quietly shrink the
universe of every enrichment computed afterwards.

## Usage

``` r
.validateLinkTable(
  linkTable,
  regions,
  weightSource = "supplied",
  verbose = TRUE
)
```

## Arguments

- linkTable:

  Data.frame or `GRanges` with the links.

- regions:

  `GRanges` with the regions of the object.

- weightSource:

  String naming where the weights come from, recorded alongside them.
  Default: `"supplied"`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A data.frame with the validated links, carrying the `region` column used
for the matching.

## Author

Sebastian Gregoricchio
