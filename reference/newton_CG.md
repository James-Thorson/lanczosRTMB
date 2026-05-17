# Nonlinear minimizer using line search with approximated Newton solution

Nonlinear minimizer designed to use cheap Hessian-vector products using
[make_Hq](https://james-thorson.github.io/lanczosRTMB/reference/make_Hq.md).
The minimizer approximates a Newton update using truncated conjugate
gradient, while adapting a regularization parameter designed and using a
line-search for each step.

## Usage

``` r
newton_CG(
  par,
  fn,
  gr,
  Hq,
  gr_tol = 1e-08,
  e_ratio = 0.01,
  maxit_newton = 100,
  maxit_CG = min(100, length(par)),
  c1 = 0.01,
  beta = 0.5,
  line_steps = 20,
  smartsearch = TRUE,
  ustep = 1,
  power = 0.5,
  u0 = 1e-04,
  stop_if_nonPD = smartsearch,
  tol10 = 0.001,
  diagnostics = FALSE,
  silent = FALSE
)
```

## Arguments

- par:

  initial parameter vector

- fn:

  function to evaluate negative log-likelihood

- gr:

  function to evaluate gradient of negative log-likelihood

- Hq:

  efficient Hessian-vector product function, e.g.,
  [make_Hq](https://james-thorson.github.io/lanczosRTMB/reference/make_Hq.md)

- gr_tol:

  early stopping condition for gradient of Newton solver

- e_ratio:

  early stopping condition for error of CG, relative to initial gradient

- maxit_newton:

  maximum iterations for Newton solver

- maxit_CG:

  maximum iterations for CG solution for each Newton iteration

- c1:

  stopping condition for line search given CG solution in each Newton
  iteration, where \\1\>c1\>0\\ seeks to prevent overshoot and resulting
  zigzagging.

- beta:

  updates in line search stepsize alpha when Armijo sufficient decrease
  condition fails

- line_steps:

  number of steps to explore for linear-search given Newton update

- smartsearch:

  Turn on adaptive stepsize algorithm for non-convex problems?

- ustep:

  Adaptive stepsize initial guess between 0 and 1.

- power:

  Parameter controlling adaptive stepsize.

- u0:

  Parameter controlling adaptive stepsize.

- stop_if_nonPD:

  whether to stop CG if recusion is not positive definite

- tol10:

  Try to exit if last 10 iterations not improved more than this.

- diagnostics:

  whether to provide extra diagnostics for each Newton iteration

- silent:

  Be silent or print progress?

## Details

This minimizer approximates Newton steps \\x\_{i+1} = x\_{i} -
H(x_i)^{-1} g(x_i)\\, where \\H(x_i) = \nabla^2 f(x_i)\\ is the Hessian
matrix and \\g(x_i) = \nabla f(x_i)\\ is the gradient of the negative
log-likelihood \\f(x_i)\\. It then uses several strategies for numerical
and computational efficiency.

1.  It approximates each Newton step \\H(x_i)^{-1} g(x_i)\\ using a
    truncated conjugate gradient algorithm that involves Hessian-vector
    products \\H(x_i) v\\. It then recursively improves the solution
    until a desired accuracy is reached, controlled by `gr_tol`, without
    ever constructing the Hessian matrix itself;

2.  It approximates each Newton solution \\H(x_i)^{-1} g(x_i)\\ without
    ever constructing or storing the Hessian \\H(x_i)\\, and instead
    computing Hessian-vector product \\H(x_i) v = \nabla( v^T \nabla
    f(x_i) )\\ using autodifferentiation in RTMB;

3.  After approximating each step direction \\H(x_i)^{-1} g(x_i)\\, it
    uses a line-search algorithm while decreasing the step-size if
    needed to find a suitable decrease in the objective function,
    controlled `line_steps`, `beta`, and `c1`;

4.  The CG is monitored to identify whether the Hessian matrix is
    positive definite (PD), and if it is not PD then the CG may be
    terminated, controlled by `stop_if_nonPD`.

5.  When using `smartsearch = TRUE`, if the Hessian is not PD, then a
    regularization `ustep` is decreased, corresponding to an increase in
    \\t=\frac{1}{u}-1\\ where the Newton step is actually using \\(H +
    tI)^{-1} g\\ where \\g\\ is the gradient, and \\1\>u\>0\\
    corresponds to \\t\>0\\. This increase or decrease in `ustep` (and
    associated decrease/increase in \\t\\) is copied from
    [`TMB::newton()`](https://rdrr.io/pkg/TMB/man/newton.html). If
    `smartsearch = FALSE` then \\ustep=1\\ and \\t=0\\ such that CG uses
    the Hessian corresponding to unregularized Newton steps. This
    smartsearch behavior controlled by `u0`, `ustep`, `power`, and
    `tol10`, and it corresponds to an adaptive "trust-region".

## Examples

``` r
library(Matrix)
library(RTMB)

# Settings
n = 10^3
rho = 0.99

# Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
P = bandSparse( n = n, k = c(-1), diagonals = list(rep(1,n)) )
Q = (Diagonal(n) - rho * t(P)) %*% (Diagonal(n) - rho * P)
x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
y = rpois( n = n, lambda = exp(1 + x) )
which_seen = sample( seq_len(n), size = n/10, replace = FALSE)
y[-which_seen] = NA

nll = function(p){
  -dgmrf(p$x, Q = Q, log = TRUE) - sum(dpois(y, lambda = exp(1 + p$x), log=TRUE), na.rm=TRUE)
}
parlist = list( x=rnorm(n) )

tape = MakeTape(nll, parlist)
gr = tape$jacfun()
Hq = make_Hq( tape, x = unlist(parlist) )
H = gr$jacfun(sparse = TRUE)

start_time = Sys.time()
opt1 = optim(
  par = unlist(parlist),
  fn = tape,
  gr = gr,
  method = "L-BFGS-B",
  control = list(
    maxit = 1e4
  )
)
opt1$runtime = Sys.time() - start_time

opt2 = newton_CG(
  par = unlist(parlist),
  fn = tape,
  gr = gr,
  Hq = Hq,
  maxit_newton = 1e4,
  gr_tol = 1e-4
)
#> value: 69002495 mgc: 5423092 ustep: 1 
#> value: 62088129 mgc: 13710652 ustep: 1 
#> value: 59887414 mgc: 47119687 ustep: 1 
#> value: 23615103 mgc: 16047732 ustep: 1 
#> value: 13011689 mgc: 12388215 ustep: 1 
#> value: 5185825 mgc: 3472569 ustep: 1 
#> value: 3855679 mgc: 3258158 ustep: 1 
#> value: 1561687 mgc: 669144.7 ustep: 1 
#> value: 886937.5 mgc: 663927.3 ustep: 1 
#> value: 592567.2 mgc: 151003.6 ustep: 1 
#> value: 360201.1 mgc: 150414 ustep: 1 
#> value: 154203.4 mgc: 141068 ustep: 1 
#> value: 76910.59 mgc: 73580.35 ustep: 1 
#> value: 49630.06 mgc: 7802.055 ustep: 1 
#> value: 24883.5 mgc: 7316.349 ustep: 1 
#> value: 22031.42 mgc: 7907.618 ustep: 1 
#> value: 9340.313 mgc: 2979.117 ustep: 1 
#> value: 5899.55 mgc: 1764.597 ustep: 1 
#> value: 3992.803 mgc: 502.0657 ustep: 1 
#> value: 2777.583 mgc: 660.3848 ustep: 1 
#> value: 2308.115 mgc: 322.1401 ustep: 1 
#> value: 1688.901 mgc: 138.9851 ustep: 1 
#> value: 1525.79 mgc: 53.88632 ustep: 1 
#> value: 1471.817 mgc: 7.037167 ustep: 1 
#> value: 1445.501 mgc: 44.51218 ustep: 1 
#> value: 1425.355 mgc: 716.6695 ustep: 1 
#> value: 1425.291 mgc: 18.56371 ustep: 1 
#> value: 1409.477 mgc: 28.93731 ustep: 1 
#> value: 1396.167 mgc: 53.65045 ustep: 1 
#> value: 1384.905 mgc: 229.0825 ustep: 1 
#> value: 1384.831 mgc: 7.669506 ustep: 1 
#> value: 1375.103 mgc: 79.23473 ustep: 1 
#> value: 1367.607 mgc: 3.710414 ustep: 1 
#> value: 1358.486 mgc: 98.72886 ustep: 1 
#> value: 1350.68 mgc: 64.2792 ustep: 1 
#> value: 1343.844 mgc: 9.320342 ustep: 1 
#> value: 1337.793 mgc: 8.665831 ustep: 1 
#> value: 1332.897 mgc: 230.3057 ustep: 1 
#> value: 1332.825 mgc: 3.615517 ustep: 1 
#> value: 1328.162 mgc: 62.29262 ustep: 1 
#> value: 1327.772 mgc: 2.776267 ustep: 1 
#> value: 1322.744 mgc: 105.3076 ustep: 1 
#> value: 1321.63 mgc: 6.371372 ustep: 1 
#> value: 1317.912 mgc: 16.60089 ustep: 1 
#> value: 1314.588 mgc: 453.7064 ustep: 1 
#> value: 1314.569 mgc: 12.83325 ustep: 1 
#> value: 1311.52 mgc: 13.30162 ustep: 1 
#> value: 1308.639 mgc: 11.31073 ustep: 1 
#> value: 1306.003 mgc: 11.31955 ustep: 1 
#> value: 1303.707 mgc: 7.989277 ustep: 1 
#> value: 1301.456 mgc: 12.12161 ustep: 1 
#> value: 1299.436 mgc: 36.06842 ustep: 1 
#> value: 1297.674 mgc: 51.79927 ustep: 1 
#> value: 1297.447 mgc: 1.616824 ustep: 1 
#> value: 1295.422 mgc: 7.32246 ustep: 1 
#> value: 1293.578 mgc: 14.8865 ustep: 1 
#> value: 1291.942 mgc: 11.15654 ustep: 1 
#> value: 1290.442 mgc: 12.67298 ustep: 1 
#> value: 1289.232 mgc: 5.237035 ustep: 1 
#> value: 1288.108 mgc: 1.313293 ustep: 1 
#> value: 1287.01 mgc: 25.87518 ustep: 1 
#> value: 1285.993 mgc: 2.945176 ustep: 1 
#> value: 1285.109 mgc: 13.23904 ustep: 1 
#> value: 1284.258 mgc: 8.969918 ustep: 1 
#> value: 1283.46 mgc: 144.8279 ustep: 1 
#> value: 1283.457 mgc: 2.171368 ustep: 1 
#> value: 1282.712 mgc: 5.880579 ustep: 1 
#> value: 1282.018 mgc: 18.63309 ustep: 1 
#> value: 1281.333 mgc: 16.00913 ustep: 1 
#> value: 1280.701 mgc: 3.761866 ustep: 1 
#> value: 1280.124 mgc: 7.287855 ustep: 1 
#> value: 1279.591 mgc: 18.0435 ustep: 1 
#> value: 1279.054 mgc: 6.61677 ustep: 1 
#> value: 1278.554 mgc: 15.57544 ustep: 1 
#> value: 1278.091 mgc: 9.883806 ustep: 1 
#> value: 1277.652 mgc: 27.47137 ustep: 1 
#> value: 1277.648 mgc: 1.387491 ustep: 1 
#> value: 1277.237 mgc: 5.397217 ustep: 1 
#> value: 1276.847 mgc: 24.92788 ustep: 1 
#> value: 1276.545 mgc: 1.323897 ustep: 1 
#> value: 1271.986 mgc: 33.73775 ustep: 1 
#> value: 1270.26 mgc: 1.999734 ustep: 1 
#> value: 1270.144 mgc: 13.15872 ustep: 1 
#> value: 1270.14 mgc: 0.7204035 ustep: 1 
#> value: 1270.045 mgc: 3.000105 ustep: 1 
#> value: 1269.955 mgc: 13.64413 ustep: 1 
#> value: 1269.955 mgc: 0.5588656 ustep: 1 
#> value: 1269.879 mgc: 12.71953 ustep: 1 
#> value: 1269.859 mgc: 0.4830983 ustep: 1 
#> value: 1269.744 mgc: 80.53164 ustep: 1 
#> value: 1269.743 mgc: 1.947211 ustep: 1 
#> value: 1269.636 mgc: 0.4455192 ustep: 1 
#> value: 1269.532 mgc: 48.30398 ustep: 1 
#> value: 1269.53 mgc: 1.530554 ustep: 1 
#> value: 1269.435 mgc: 4.408405 ustep: 1 
#> value: 1269.336 mgc: 3.524596 ustep: 1 
#> value: 1269.244 mgc: 3.964306 ustep: 1 
#> value: 1269.138 mgc: 21.49671 ustep: 1 
#> value: 1269.128 mgc: 1.034898 ustep: 1 
#> value: 1269.04 mgc: 4.932344 ustep: 1 
#> value: 1269.038 mgc: 0.1879522 ustep: 1 
#> value: 1268.781 mgc: 13.78388 ustep: 1 
#> value: 1268.607 mgc: 0.8582689 ustep: 1 
#> value: 1268.574 mgc: 3.428539 ustep: 1 
#> value: 1268.547 mgc: 10.65905 ustep: 1 
#> value: 1268.547 mgc: 0.4780089 ustep: 1 
#> value: 1268.521 mgc: 2.412192 ustep: 1 
#> value: 1268.497 mgc: 0.3939968 ustep: 1 
#> value: 1268.473 mgc: 8.043495 ustep: 1 
#> value: 1268.47 mgc: 0.2732363 ustep: 1 
#> value: 1268.44 mgc: 18.37112 ustep: 1 
#> value: 1268.434 mgc: 1.142214 ustep: 1 
#> value: 1268.412 mgc: 2.029332 ustep: 1 
#> value: 1268.388 mgc: 3.180423 ustep: 1 
#> value: 1268.365 mgc: 1.242234 ustep: 1 
#> value: 1268.345 mgc: 0.5126696 ustep: 1 
#> value: 1268.324 mgc: 0.4391103 ustep: 1 
#> value: 1268.305 mgc: 4.543721 ustep: 1 
#> value: 1268.292 mgc: 0.230684 ustep: 1 
#> value: 1268.228 mgc: 50.42631 ustep: 1 
#> value: 1268.228 mgc: 2.551755 ustep: 1 
#> value: 1268.172 mgc: 0.3338634 ustep: 1 
#> value: 1268.126 mgc: 6.026182 ustep: 1 
#> value: 1268.088 mgc: 0.2713439 ustep: 1 
#> value: 1268.035 mgc: 10.39659 ustep: 1 
#> value: 1268.025 mgc: 0.4893303 ustep: 1 
#> value: 1267.988 mgc: 1.080293 ustep: 1 
#> value: 1267.959 mgc: 1.919685 ustep: 1 
#> value: 1267.932 mgc: 2.455357 ustep: 1 
#> value: 1267.911 mgc: 10.5314 ustep: 1 
#> value: 1267.903 mgc: 0.5010142 ustep: 1 
#> value: 1267.89 mgc: 11.15183 ustep: 1 
#> value: 1267.889 mgc: 0.4864803 ustep: 1 
#> value: 1267.878 mgc: 1.825475 ustep: 1 
#> value: 1267.867 mgc: 17.76147 ustep: 1 
#> value: 1267.866 mgc: 1.645518 ustep: 1 
#> value: 1267.857 mgc: 1.885765 ustep: 1 
#> value: 1267.848 mgc: 0.1500001 ustep: 1 
#> value: 1267.836 mgc: 0.9098671 ustep: 1 
#> value: 1267.824 mgc: 0.190068 ustep: 1 
#> value: 1267.813 mgc: 9.423086 ustep: 1 
#> value: 1267.812 mgc: 0.458711 ustep: 1 
#> value: 1267.802 mgc: 0.9599396 ustep: 1 
#> value: 1267.792 mgc: 22.61884 ustep: 1 
#> value: 1267.792 mgc: 0.2377696 ustep: 1 
#> value: 1267.783 mgc: 0.6577689 ustep: 1 
#> value: 1267.774 mgc: 0.33732 ustep: 1 
#> value: 1267.765 mgc: 4.81276 ustep: 1 
#> value: 1267.765 mgc: 0.05260067 ustep: 1 
#> value: 1267.754 mgc: 0.4577919 ustep: 1 
#> value: 1267.745 mgc: 2.220954 ustep: 1 
#> value: 1267.745 mgc: 0.0641027 ustep: 1 
#> value: 1267.734 mgc: 1.946944 ustep: 1 
#> value: 1267.724 mgc: 0.4196361 ustep: 1 
#> value: 1267.716 mgc: 6.087933 ustep: 1 
#> value: 1267.716 mgc: 0.6407287 ustep: 1 
#> value: 1267.709 mgc: 13.84851 ustep: 1 
#> value: 1267.709 mgc: 0.7590716 ustep: 1 
#> value: 1267.701 mgc: 3.325148 ustep: 1 
#> value: 1267.694 mgc: 0.1835423 ustep: 1 
#> value: 1267.684 mgc: 0.7930306 ustep: 1 
#> value: 1267.674 mgc: 0.4343502 ustep: 1 
#> value: 1267.664 mgc: 33.18544 ustep: 1 
#> value: 1267.664 mgc: 0.5055578 ustep: 1 
#> value: 1267.655 mgc: 1.541604 ustep: 1 
#> value: 1267.647 mgc: 0.5222939 ustep: 1 
#> value: 1267.64 mgc: 0.3163822 ustep: 1 
#> value: 1267.634 mgc: 1.425528 ustep: 1 
#> value: 1267.628 mgc: 0.1411388 ustep: 1 
#> value: 1267.622 mgc: 0.4182761 ustep: 1 
#> value: 1267.617 mgc: 0.5258442 ustep: 1 
#> value: 1267.613 mgc: 2.854774 ustep: 1 
#> value: 1267.609 mgc: 0.1980784 ustep: 1 
#> value: 1267.597 mgc: 0.8861654 ustep: 1 
#> value: 1267.587 mgc: 2.807737 ustep: 1 
#> value: 1267.586 mgc: 0.1681381 ustep: 1 
#> value: 1267.583 mgc: 3.569363 ustep: 1 
#> value: 1267.583 mgc: 0.1180857 ustep: 1 
#> value: 1267.578 mgc: 0.321002 ustep: 1 
#> value: 1267.573 mgc: 0.5474265 ustep: 1 
#> value: 1267.569 mgc: 43.35063 ustep: 1 
#> value: 1267.569 mgc: 0.1704768 ustep: 1 
#> value: 1267.566 mgc: 0.9920958 ustep: 1 
#> value: 1267.563 mgc: 3.101464 ustep: 1 
#> value: 1267.563 mgc: 0.2371388 ustep: 1 
#> value: 1267.561 mgc: 0.5493272 ustep: 1 
#> value: 1267.558 mgc: 0.6097116 ustep: 1 
#> value: 1267.556 mgc: 0.4706305 ustep: 1 
#> value: 1267.554 mgc: 9.366882 ustep: 1 
#> value: 1267.554 mgc: 0.5738787 ustep: 1 
#> value: 1267.552 mgc: 0.1528771 ustep: 1 
#> value: 1267.55 mgc: 0.1335229 ustep: 1 
#> value: 1267.548 mgc: 0.5279556 ustep: 1 
#> value: 1267.546 mgc: 6.489676 ustep: 1 
#> value: 1267.546 mgc: 0.5906358 ustep: 1 
#> value: 1267.544 mgc: 6.970496 ustep: 1 
#> value: 1267.544 mgc: 0.516262 ustep: 1 
#> value: 1267.542 mgc: 0.2784254 ustep: 1 
#> value: 1267.539 mgc: 1.161649 ustep: 1 
#> value: 1267.537 mgc: 3.88863 ustep: 1 
#> value: 1267.537 mgc: 0.09139713 ustep: 1 
#> value: 1267.535 mgc: 0.9254993 ustep: 1 
#> value: 1267.533 mgc: 0.1457956 ustep: 1 
#> value: 1267.531 mgc: 0.4526357 ustep: 1 
#> value: 1267.529 mgc: 0.2739218 ustep: 1 
#> value: 1267.526 mgc: 0.1988614 ustep: 1 
#> value: 1267.524 mgc: 1.26461 ustep: 1 
#> value: 1267.523 mgc: 0.6468486 ustep: 1 
#> value: 1267.521 mgc: 0.3147832 ustep: 1 
#> value: 1267.519 mgc: 1.359603 ustep: 1 
#> value: 1267.517 mgc: 0.06082461 ustep: 1 
#> value: 1267.51 mgc: 3.944497 ustep: 1 
#> value: 1267.506 mgc: 0.1757769 ustep: 1 
#> value: 1267.503 mgc: 8.378352 ustep: 1 
#> value: 1267.503 mgc: 0.7946084 ustep: 1 
#> value: 1267.5 mgc: 0.2005704 ustep: 1 
#> value: 1267.498 mgc: 2.766794 ustep: 1 
#> value: 1267.498 mgc: 0.06094397 ustep: 1 
#> value: 1267.495 mgc: 2.690274 ustep: 1 
#> value: 1267.495 mgc: 0.1045404 ustep: 1 
#> value: 1267.493 mgc: 0.2379606 ustep: 1 
#> value: 1267.491 mgc: 0.9025408 ustep: 1 
#> value: 1267.488 mgc: 0.7355307 ustep: 1 
#> value: 1267.486 mgc: 0.1450908 ustep: 1 
#> value: 1267.484 mgc: 3.370561 ustep: 1 
#> value: 1267.484 mgc: 0.2565902 ustep: 1 
#> value: 1267.482 mgc: 1.017969 ustep: 1 
#> value: 1267.479 mgc: 0.2123067 ustep: 1 
#> value: 1267.477 mgc: 0.7938712 ustep: 1 
#> value: 1267.475 mgc: 0.3612189 ustep: 1 
#> value: 1267.473 mgc: 1.293369 ustep: 1 
#> value: 1267.471 mgc: 1.875097 ustep: 1 
#> value: 1267.471 mgc: 0.04682026 ustep: 1 
#> value: 1267.469 mgc: 0.5072438 ustep: 1 
#> value: 1267.467 mgc: 0.07069927 ustep: 1 
#> value: 1267.465 mgc: 1.088561 ustep: 1 
#> value: 1267.463 mgc: 0.04207041 ustep: 1 
#> value: 1267.443 mgc: 36.0518 ustep: 1 
#> value: 1267.443 mgc: 2.730045 ustep: 1 
#> value: 1267.427 mgc: 0.08506821 ustep: 1 
#> value: 1267.425 mgc: 1.3402 ustep: 1 
#> value: 1267.425 mgc: 0.08655548 ustep: 1 
#> value: 1267.422 mgc: 0.7728272 ustep: 1 
#> value: 1267.42 mgc: 0.04267999 ustep: 1 
#> value: 1267.418 mgc: 0.6294577 ustep: 1 
#> value: 1267.416 mgc: 0.04979805 ustep: 1 
#> value: 1267.406 mgc: 1.597226 ustep: 1 
#> value: 1267.397 mgc: 0.07300517 ustep: 1 
#> value: 1267.395 mgc: 0.9721355 ustep: 1 
#> value: 1267.394 mgc: 0.04164901 ustep: 1 
#> value: 1267.392 mgc: 2.677317 ustep: 1 
#> value: 1267.392 mgc: 0.04846821 ustep: 1 
#> value: 1267.39 mgc: 2.27303 ustep: 1 
#> value: 1267.39 mgc: 0.1135402 ustep: 1 
#> value: 1267.389 mgc: 0.9086769 ustep: 1 
#> value: 1267.388 mgc: 0.04185761 ustep: 1 
#> value: 1267.386 mgc: 2.716847 ustep: 1 
#> value: 1267.386 mgc: 0.1276475 ustep: 1 
#> value: 1267.385 mgc: 2.466907 ustep: 1 
#> value: 1267.385 mgc: 0.04846515 ustep: 1 
#> value: 1267.384 mgc: 1.274862 ustep: 1 
#> value: 1267.382 mgc: 0.08011217 ustep: 1 
#> value: 1267.381 mgc: 1.121125 ustep: 1 
#> value: 1267.38 mgc: 0.07026672 ustep: 1 
#> value: 1267.379 mgc: 0.1262073 ustep: 1 
#> value: 1267.378 mgc: 0.0659245 ustep: 1 
#> value: 1267.378 mgc: 1.379214 ustep: 1 
#> value: 1267.378 mgc: 0.06503061 ustep: 1 
#> value: 1267.377 mgc: 0.4380098 ustep: 1 
#> value: 1267.376 mgc: 0.2328043 ustep: 1 
#> value: 1267.375 mgc: 0.2463847 ustep: 1 
#> value: 1267.375 mgc: 1.95885 ustep: 1 
#> value: 1267.375 mgc: 0.05691374 ustep: 1 
#> value: 1267.374 mgc: 0.1390638 ustep: 1 
#> value: 1267.373 mgc: 0.5053454 ustep: 1 
#> value: 1267.373 mgc: 1.139701 ustep: 1 
#> value: 1267.373 mgc: 0.06574726 ustep: 1 
#> value: 1267.368 mgc: 39.28754 ustep: 1 
#> value: 1267.368 mgc: 1.310191 ustep: 1 
#> value: 1267.365 mgc: 0.05627178 ustep: 1 
#> value: 1267.364 mgc: 2.406226 ustep: 1 
#> value: 1267.363 mgc: 0.1281036 ustep: 1 
#> value: 1267.363 mgc: 0.7517985 ustep: 1 
#> value: 1267.362 mgc: 0.01913347 ustep: 1 
#> value: 1267.361 mgc: 0.2758798 ustep: 1 
#> value: 1267.36 mgc: 0.110571 ustep: 1 
#> value: 1267.359 mgc: 1.859292 ustep: 1 
#> value: 1267.359 mgc: 0.0668008 ustep: 1 
#> value: 1267.358 mgc: 0.1274225 ustep: 1 
#> value: 1267.357 mgc: 1.179774 ustep: 1 
#> value: 1267.357 mgc: 0.05772288 ustep: 1 
#> value: 1267.357 mgc: 7.103074 ustep: 1 
#> value: 1267.357 mgc: 0.2326094 ustep: 1 
#> value: 1267.356 mgc: 0.5548701 ustep: 1 
#> value: 1267.355 mgc: 0.169072 ustep: 1 
#> value: 1267.355 mgc: 1.384649 ustep: 1 
#> value: 1267.355 mgc: 0.04965044 ustep: 1 
#> value: 1267.354 mgc: 5.87694 ustep: 1 
#> value: 1267.354 mgc: 0.1123729 ustep: 1 
#> value: 1267.354 mgc: 1.695678 ustep: 1 
#> value: 1267.354 mgc: 0.0464295 ustep: 1 
#> value: 1267.353 mgc: 1.505321 ustep: 1 
#> value: 1267.353 mgc: 0.04147998 ustep: 1 
#> value: 1267.352 mgc: 0.8458953 ustep: 1 
#> value: 1267.352 mgc: 0.02957149 ustep: 1 
#> value: 1267.351 mgc: 5.420434 ustep: 1 
#> value: 1267.351 mgc: 0.3100229 ustep: 1 
#> value: 1267.351 mgc: 0.4006729 ustep: 1 
#> value: 1267.35 mgc: 0.4740931 ustep: 1 
#> value: 1267.349 mgc: 0.05979546 ustep: 1 
#> value: 1267.349 mgc: 0.5049684 ustep: 1 
#> value: 1267.348 mgc: 0.7675929 ustep: 1 
#> value: 1267.347 mgc: 3.224137 ustep: 1 
#> value: 1267.347 mgc: 0.2078309 ustep: 1 
#> value: 1267.347 mgc: 0.1154254 ustep: 1 
#> value: 1267.346 mgc: 0.1853998 ustep: 1 
#> value: 1267.346 mgc: 0.5521956 ustep: 1 
#> value: 1267.345 mgc: 0.8802601 ustep: 1 
#> value: 1267.344 mgc: 0.05039581 ustep: 1 
#> value: 1267.344 mgc: 0.5070966 ustep: 1 
#> value: 1267.343 mgc: 0.03091531 ustep: 1 
#> value: 1267.343 mgc: 0.2907092 ustep: 1 
#> value: 1267.342 mgc: 0.7324482 ustep: 1 
#> value: 1267.341 mgc: 36.00393 ustep: 1 
#> value: 1267.341 mgc: 2.576899 ustep: 1 
#> value: 1267.341 mgc: 0.1818308 ustep: 1 
#> value: 1267.34 mgc: 0.03196487 ustep: 1 
#> value: 1267.34 mgc: 0.4306464 ustep: 1 
#> value: 1267.339 mgc: 1.678724 ustep: 1 
#> value: 1267.339 mgc: 0.07937939 ustep: 1 
#> value: 1267.339 mgc: 1.965848 ustep: 1 
#> value: 1267.339 mgc: 0.07790287 ustep: 1 
#> value: 1267.338 mgc: 0.1415751 ustep: 1 
#> value: 1267.338 mgc: 15.18634 ustep: 1 
#> value: 1267.338 mgc: 0.2223445 ustep: 1 
#> value: 1267.337 mgc: 0.178543 ustep: 1 
#> value: 1267.337 mgc: 3.361097 ustep: 1 
#> value: 1267.337 mgc: 0.0543379 ustep: 1 
#> value: 1267.336 mgc: 0.3863202 ustep: 1 
#> value: 1267.336 mgc: 0.2066731 ustep: 1 
#> value: 1267.335 mgc: 2.723113 ustep: 1 
#> value: 1267.335 mgc: 0.06216811 ustep: 1 
#> value: 1267.335 mgc: 1.154198 ustep: 1 
#> value: 1267.335 mgc: 0.05110498 ustep: 1 
#> value: 1267.334 mgc: 0.1305264 ustep: 1 
#> value: 1267.334 mgc: 0.4646031 ustep: 1 
#> value: 1267.333 mgc: 0.1013669 ustep: 1 
#> value: 1267.333 mgc: 0.3248838 ustep: 1 
#> value: 1267.332 mgc: 1.332741 ustep: 1 
#> value: 1267.332 mgc: 0.0792638 ustep: 1 
#> value: 1267.332 mgc: 0.2469052 ustep: 1 
#> value: 1267.332 mgc: 0.08218199 ustep: 1 
#> value: 1267.331 mgc: 0.2763074 ustep: 1 
#> value: 1267.331 mgc: 2.456232 ustep: 1 
#> value: 1267.331 mgc: 0.1203033 ustep: 1 
#> value: 1267.33 mgc: 0.4327469 ustep: 1 
#> value: 1267.33 mgc: 0.2734012 ustep: 1 
#> value: 1267.33 mgc: 2.065591 ustep: 1 
#> value: 1267.33 mgc: 0.1623543 ustep: 1 
#> value: 1267.329 mgc: 0.2492188 ustep: 1 
#> value: 1267.329 mgc: 1.022884 ustep: 1 
#> value: 1267.329 mgc: 0.06033905 ustep: 1 
#> value: 1267.328 mgc: 4.735754 ustep: 1 
#> value: 1267.328 mgc: 0.1191048 ustep: 1 
#> value: 1267.328 mgc: 0.2145898 ustep: 1 
#> value: 1267.327 mgc: 3.315181 ustep: 1 
#> value: 1267.327 mgc: 0.1034758 ustep: 1 
#> value: 1267.327 mgc: 0.09816467 ustep: 1 
#> value: 1267.326 mgc: 0.1833721 ustep: 1 
#> value: 1267.326 mgc: 0.09216511 ustep: 1 
#> value: 1267.326 mgc: 0.05990282 ustep: 1 
#> value: 1267.325 mgc: 0.7880752 ustep: 1 
#> value: 1267.325 mgc: 0.03315138 ustep: 1 
#> value: 1267.325 mgc: 0.2601117 ustep: 1 
#> value: 1267.324 mgc: 0.3478648 ustep: 1 
#> value: 1267.324 mgc: 5.224672 ustep: 1 
#> value: 1267.324 mgc: 0.3186263 ustep: 1 
#> value: 1267.324 mgc: 0.05951982 ustep: 1 
#> value: 1267.323 mgc: 0.1619244 ustep: 1 
#> value: 1267.323 mgc: 0.1021333 ustep: 1 
#> value: 1267.322 mgc: 1.75757 ustep: 1 
#> value: 1267.322 mgc: 0.0589379 ustep: 1 
#> value: 1267.321 mgc: 0.1543897 ustep: 1 
#> value: 1267.321 mgc: 1.595594 ustep: 1 
#> value: 1267.321 mgc: 0.04362338 ustep: 1 
#> value: 1267.321 mgc: 0.09973457 ustep: 1 
#> value: 1267.32 mgc: 3.320548 ustep: 1 
#> value: 1267.32 mgc: 0.1465666 ustep: 1 
#> value: 1267.32 mgc: 0.275159 ustep: 1 
#> value: 1267.32 mgc: 2.444444 ustep: 1 
#> value: 1267.319 mgc: 0.169374 ustep: 1 
#> value: 1267.319 mgc: 0.8856646 ustep: 1 
#> value: 1267.319 mgc: 0.03836001 ustep: 1 
#> value: 1267.319 mgc: 1.415484 ustep: 1 
#> value: 1267.318 mgc: 0.05771761 ustep: 1 
#> value: 1267.318 mgc: 0.1841386 ustep: 1 
#> value: 1267.318 mgc: 3.577176 ustep: 1 
#> value: 1267.318 mgc: 0.06269276 ustep: 1 
#> value: 1267.317 mgc: 0.8460506 ustep: 1 
#> value: 1267.317 mgc: 0.05341851 ustep: 1 
#> value: 1267.317 mgc: 1.017335 ustep: 1 
#> value: 1267.317 mgc: 0.03925536 ustep: 1 
#> value: 1267.316 mgc: 0.1439491 ustep: 1 
#> value: 1267.316 mgc: 0.1193836 ustep: 1 
#> value: 1267.316 mgc: 0.06187685 ustep: 1 
#> value: 1267.315 mgc: 0.5210019 ustep: 1 
#> value: 1267.315 mgc: 0.05018722 ustep: 1 
#> value: 1267.315 mgc: 0.0961183 ustep: 1 
#> value: 1267.314 mgc: 0.205766 ustep: 1 
#> value: 1267.314 mgc: 3.55181 ustep: 1 
#> value: 1267.314 mgc: 0.1620508 ustep: 1 
#> value: 1267.314 mgc: 0.2633271 ustep: 1 
#> value: 1267.313 mgc: 0.6460342 ustep: 1 
#> value: 1267.313 mgc: 0.04082085 ustep: 1 
#> value: 1267.313 mgc: 2.848909 ustep: 1 
#> value: 1267.313 mgc: 0.05326432 ustep: 1 
#> value: 1267.312 mgc: 0.06965334 ustep: 1 
#> value: 1267.312 mgc: 0.07558166 ustep: 1 
#> value: 1267.312 mgc: 0.6273978 ustep: 1 
#> value: 1267.312 mgc: 0.02927192 ustep: 1 
#> value: 1267.312 mgc: 0.1776044 ustep: 1 
#> value: 1267.311 mgc: 0.2524507 ustep: 1 
#> value: 1267.311 mgc: 0.5714615 ustep: 1 
#> value: 1267.311 mgc: 0.03634663 ustep: 1 
#> value: 1267.31 mgc: 6.569855 ustep: 1 
#> value: 1267.31 mgc: 0.1836578 ustep: 1 
#> value: 1267.31 mgc: 0.4845608 ustep: 1 
#> value: 1267.31 mgc: 0.2200542 ustep: 1 
#> value: 1267.309 mgc: 0.08006179 ustep: 1 
#> value: 1267.309 mgc: 2.06136 ustep: 1 
#> value: 1267.309 mgc: 0.06609773 ustep: 1 
#> value: 1267.309 mgc: 0.1448974 ustep: 1 
#> value: 1267.309 mgc: 0.06329128 ustep: 1 
#> value: 1267.308 mgc: 1.120575 ustep: 1 
#> value: 1267.308 mgc: 0.055573 ustep: 1 
#> value: 1267.308 mgc: 0.5240807 ustep: 1 
#> value: 1267.308 mgc: 0.8907478 ustep: 1 
#> value: 1267.308 mgc: 0.04499826 ustep: 1 
#> value: 1267.308 mgc: 0.503043 ustep: 1 
#> value: 1267.307 mgc: 0.02490453 ustep: 1 
#> value: 1267.307 mgc: 0.09966598 ustep: 1 
#> value: 1267.306 mgc: 0.02953167 ustep: 1 
#> value: 1267.306 mgc: 0.1516249 ustep: 1 
#> value: 1267.306 mgc: 0.03704595 ustep: 1 
#> value: 1267.305 mgc: 1.26063 ustep: 1 
#> value: 1267.305 mgc: 0.04350208 ustep: 1 
#> value: 1267.305 mgc: 0.4473574 ustep: 1 
#> value: 1267.304 mgc: 1.945948 ustep: 1 
#> value: 1267.304 mgc: 0.06838389 ustep: 1 
#> value: 1267.304 mgc: 0.8137347 ustep: 1 
#> value: 1267.304 mgc: 0.04319604 ustep: 1 
#> value: 1267.304 mgc: 0.3106419 ustep: 1 
#> value: 1267.304 mgc: 0.1588633 ustep: 1 
#> value: 1267.303 mgc: 0.1902695 ustep: 1 
#> value: 1267.303 mgc: 0.2218362 ustep: 1 
#> value: 1267.303 mgc: 0.1798331 ustep: 1 
#> value: 1267.303 mgc: 0.1387916 ustep: 1 
#> value: 1267.302 mgc: 17.63813 ustep: 1 
#> value: 1267.302 mgc: 0.6930635 ustep: 1 
#> value: 1267.302 mgc: 0.03094613 ustep: 1 
#> value: 1267.299 mgc: 6.404396 ustep: 1 
#> value: 1267.298 mgc: 0.4322432 ustep: 1 
#> value: 1267.296 mgc: 0.03137295 ustep: 1 
#> value: 1267.294 mgc: 0.6520718 ustep: 1 
#> value: 1267.292 mgc: 0.04491016 ustep: 1 
#> value: 1267.29 mgc: 3.008786 ustep: 1 
#> value: 1267.29 mgc: 0.1438102 ustep: 1 
#> value: 1267.288 mgc: 0.0311361 ustep: 1 
#> value: 1267.287 mgc: 0.1132093 ustep: 1 
#> value: 1267.286 mgc: 0.1403634 ustep: 1 
#> value: 1267.286 mgc: 0.03119828 ustep: 1 
#> value: 1267.286 mgc: 0.8872719 ustep: 1 
#> value: 1267.286 mgc: 0.04546892 ustep: 1 
#> value: 1267.286 mgc: 0.04113355 ustep: 1 
#> value: 1267.285 mgc: 0.1791689 ustep: 1 
#> value: 1267.285 mgc: 1.114783 ustep: 1 
#> value: 1267.285 mgc: 0.04825858 ustep: 1 
#> value: 1267.285 mgc: 0.2973368 ustep: 1 
#> value: 1267.285 mgc: 0.04309803 ustep: 1 
#> value: 1267.285 mgc: 0.1115932 ustep: 1 
#> value: 1267.284 mgc: 0.1243357 ustep: 1 
#> value: 1267.284 mgc: 4.801192 ustep: 1 
#> value: 1267.284 mgc: 0.4416367 ustep: 1 
#> value: 1267.284 mgc: 0.6995811 ustep: 1 
#> value: 1267.284 mgc: 0.03277868 ustep: 1 
#> value: 1267.284 mgc: 6.548815 ustep: 1 
#> value: 1267.284 mgc: 0.470584 ustep: 1 
#> value: 1267.284 mgc: 0.0219272 ustep: 1 
#> Not improving much - will try early exit...PD hess?: TRUE 

# Compare the estimates and speed
matplot( cbind(opt1$par, opt2$par), type = "l", col = c("black","blue","red"), lty = "solid" )

c(opt1$runtime, opt2$runtime)
#> Time differences in secs
#> [1]  5.800411 34.661608
```
