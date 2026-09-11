# RegionSetDE.setScores class

S4 class storing the signal scores computed by
[`scoreRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreRegionSets.md),
one per region set per library, together with the paired comparisons
between the sets. No contrast lies behind these numbers and no model was
fitted on them, which is what separates this class from
`RegionSetDE.setResults`: the replication is the libraries, and the
comparison runs between sets rather than between conditions.

## Slots

- `scores`:

  Data.frame with one row per library per set, carrying the summarised
  signal, the reference it was divided by and the `score` on a log2
  scale.

- `comparisons`:

  Data.frame with one row per pair of sets, carrying the mean paired
  difference between their scores, its confidence interval, the paired
  t-statistic and the adjusted p-value.

- `sample.metadata`:

  Data.frame with the `colData` of the counts object the scores came
  from, so that the libraries can be grouped and coloured downstream
  without handing that object over again.

- `reference`:

  String naming what the signal of a set was divided by, one among
  `"background"`, `"regions"` and `"none"`.

- `assay`:

  String with the name of the assay the signal was read from.

- `summary`:

  String with how the regions of a set were summarised, either `"mean"`
  or `"median"`.

- `per.basepair`:

  Logical value indicating whether the signal was divided by the width
  of the region before being summarised.

## See also

[`scoreRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreRegionSets.md)

## Author

Sebastian Gregoricchio
