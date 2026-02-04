# Check if a Processing Step has been Applied

Scans an object's processing history to determine if a step with a
specific name has already been performed. This is useful for preventing
redundant operations or for enforcing a required order of preprocessing
steps.

## Usage

``` r
has_step(object, step_name)
```

## Arguments

- object:

  An S4 object containing a `processing` slot (list).

- step_name:

  A character string specifying the unique name of the step to check for
  (e.g., `"log2_transform"`).

## Value

A logical value: `TRUE` if a log entry with the given `step_name` exists
in the object's history, `FALSE` otherwise.
