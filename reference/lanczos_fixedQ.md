# Calculate Lanczos approximation using fixed Q

Recalculate Lanczos approximation using fixed Q to speed up gradient
calculations

## Usage

``` r
lanczos_fixedQ(
  Hq,
  Q,
  x = attr(Hq, "env")$x0,
  V = matrix(1, nrow = length(x), ncol = 0)
)
```

## Arguments

- Hq:

  function that calculates the product `H %*% q` given probe `q` and
  parameters `x`

- Q:

  matrix of fixed probes

- x:

  parameter vector used when calculating the Hessian matrix

- V:

  deflation matrix, where we transform probes \\q' = (I-V V^T) q\\ and
  the Hessian by \\(I-V V^T) H (I-V V^T)\\ to eliminate axes in \\V\\
