# RegionSetDE.geneResults class

S4 class holding the per-gene table built by
[`aggregateToGenes`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/aggregateToGenes.md)
from a per-region result. One row per gene, with the regions that
reached it counted and their directions reported rather than resolved,
and the filters applied upstream still readable in the inherited
provenance slots.

## Slots

- `results`:

  Data.frame with one row per gene. Read it with
  [`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md).

- `contrast`:

  String with the contrast the regions were tested on.

- `method`:

  String with the aggregation used, `"simes"` or `"potential"`.

- `links`:

  Data.frame with the region to gene links the aggregation was built on,
  one row per link.

- `weighting`:

  String saying what weighted the regions of a gene, one among
  `"link.weight"`, `"distance"` and `"none"`.

- `engine`:

  String with the engine the regions were fitted with.

- `thresholds`:

  List with the cut-offs the directions were read at.

## Author

Sebastian Gregoricchio
