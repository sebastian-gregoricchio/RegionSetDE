# RegionSetDE.provenance class

Virtual class collecting the filters and the parameters shared by all
the RegionSetDE objects, so that the origin of the regions survives
every step of the analysis. Not meant to be instantiated directly.

## Slots

- `blacklist`:

  `GRanges` with the regions removed by
  [`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md),
  `NULL` when no blacklist has been applied.

- `whitelist`:

  `GRanges` with the regions used by
  [`applyWhitelist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyWhitelist.md)
  to restrict the sets, `NULL` when no whitelist has been applied.

- `genome.assembly`:

  String with the genome assembly of the regions, `NULL` when not
  declared.

- `seqlevels.style`:

  String with the chromosome naming style shared by all the sets.

- `filtering.log`:

  Data.frame collecting the number of regions before and after each
  filtering step.

- `parameters`:

  List with the arguments used at each step.

## Author

Sebastian Gregoricchio
