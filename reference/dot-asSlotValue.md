# .asSlotValue

Brings a value to the class of the slot it goes in. A data.frame meant
for a `DataFrame` slot is converted, so that the functions of the
package can keep building their tables with dplyr.

## Usage

``` r
.asSlotValue(className, slotName, value)
```

## Arguments

- className:

  String with the class of the object.

- slotName:

  String with the name of the slot.

- value:

  The value to store.

## Value

The value, converted when the slot needs it.

## Author

Sebastian Gregoricchio
