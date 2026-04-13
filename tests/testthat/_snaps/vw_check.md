# check_numeric_param: non-numeric errors

    Code
      check_numeric_param(x)
    Condition
      Error in `check_numeric_param()`:
      ! x must be a single numeric value.

# check_numeric_param: non-scalar errors

    Code
      check_numeric_param(x)
    Condition
      Error in `check_numeric_param()`:
      ! x must be a single numeric value.

# check_numeric_param: integer constraint errors

    Code
      check_numeric_param(x, integer = TRUE)
    Condition
      Error in `check_numeric_param()`:
      ! x must be an integer.

# check_numeric_param: range constraint errors

    Code
      check_numeric_param(x, lower = 1, upper = 10)
    Condition
      Error in `check_numeric_param()`:
      ! x must be between: 1 and 10.

---

    Code
      check_numeric_param(y, lower = 1, upper = 10)
    Condition
      Error in `check_numeric_param()`:
      ! y must be between: 1 and 10.

# check_numeric_param: set membership errors

    Code
      check_numeric_param(x, set = c(2, 3, 4))
    Condition
      Error in `check_numeric_param()`:
      ! x must be one of: 2, 3, or 4.

