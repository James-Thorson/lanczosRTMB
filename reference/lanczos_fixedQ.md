# Calculate Lanczos approximation using fixed Q

Recalculate Lanczos approximation using fixed Q to speed up gradient
calculations

## Usage

``` r
lanczos_fixedQ(Hq, Q, x = attr(Hq, "env")$x0)
```

## Arguments

- Hq:

  function that calculates the product `H %*% q` given probe `q` and
  parameters `x`

- Q:

  matrix of fixed probes

- x:

  parameter vector used when calculating the Hessian matrix
