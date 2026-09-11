# RegionSetDE.setResultsList class

S4 class holding the set level outcome of several contrasts run on the
same fit. Every function that takes a `RegionSetDE.setResults` object
also takes this one, together with a `contrast` argument naming which of
them to use.

## Slots

- `results`:

  Named list of `RegionSetDE.setResults` objects.

- `contrasts`:

  Character vector with the names of the contrasts, in the order they
  were run.

## Author

Sebastian Gregoricchio
