# Estimate diag(H) using Hutchinson probes

Estimate diag(H) using Hutchinson probes

## Usage

``` r
estimate_diag_hutchinson(
  Hq,
  n_probes = 30,
  x = attr(Hq, "env")$x0,
  n = length(attr(Hq, "env")$which_random)
)
```
