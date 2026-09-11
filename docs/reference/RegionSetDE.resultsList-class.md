# RegionSetDE.resultsList class

S4 class holding the outcome of several contrasts run on the same fit.
Since all of them come from one model, the design, the offsets and the
dispersion are shared and the contrasts can be compared with each other
directly. Every function that takes a `RegionSetDE.results` object also
takes this one, together with a `contrast` argument naming which of them
to use.

## Slots

- `results`:

  Named list of `RegionSetDE.results` objects.

- `contrasts`:

  Character vector with the names of the contrasts, in the order they
  were run.

## Author

Sebastian Gregoricchio
