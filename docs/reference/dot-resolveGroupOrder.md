# .resolveGroupOrder

Decides the order the groups appear in, putting the reference level of
the contrast on the left so that it gets the blue family.

## Usage

``` r
.resolveGroupOrder(groupOrder = NULL, groupLevels, contrastGroups = NULL)
```

## Arguments

- groupOrder:

  Character vector given by the user, or `NULL`.

- groupLevels:

  Character vector with the groups present in the object.

- contrastGroups:

  Character vector of length two with the levels the contrast compares,
  the first being the one the fold change is positive for. Default:
  `NULL`.

## Value

A character vector with every group, ordered.

## Author

Sebastian Gregoricchio
