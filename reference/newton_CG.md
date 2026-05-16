# Nonlinear minimizer using line search with approximated Newton solution

Nonlinear minimizer designed for cheap Hessian-vector products using
[make_Hq](https://james-thorson.github.io/lanczosRTMB/reference/make_Hq.md),
involving iterating a linear search along the truncated conjugate
gradient for a Newton.

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
  smartsearch = FALSE,
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
  iteration

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

This minimizer approximates Newton steps, but solving each Newton
solution using a truncated conjugate gradient using Hessian-vector
products without ever constructing the Hessian matrix directly. It then
uses three separate tricks for numerical and computational efficiency.

1.  For each CG, it recursively improves the solution until a desired
    accuracy is reached;

2.  For each Newton-step, it uses a line-search algorithm, while
    decreasing the step-size to find a decrease in objective function.

3.  The CG is monitored to identify whether the Hessian matrix is
    positive definite (PD), and if it is not PD then the CG may be
    terminated.

4.  When using `smartsearch = TRUE`, if the Hessian is not PD, then a
    regularization `ustep` is decreased, corresponding to an increase in
    `t` where the Newton step is actually using \\(H + tI)^{-1} g\\
    where \\g\\ is the gradient. This increase or decrease in `ustep`
    (and associated decrease/increase in \\t\\) is copied from
    [`TMB::newton`](https://rdrr.io/pkg/TMB/man/newton.html). If
    `smartsearch = FALSE` then \\t=0\\ and CG uses the Hessian
    corresponding to unregularized Newton steps.

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
#> value: 69002495 mgc: 5423092 
#> value: 62088129 mgc: 13710652 
#> value: 59887414 mgc: 47119687 
#> value: 23615103 mgc: 16047732 
#> value: 13011689 mgc: 12388215 
#> value: 5185825 mgc: 3472569 
#> value: 3855679 mgc: 3258158 
#> value: 1561687 mgc: 669144.7 
#> value: 886937.5 mgc: 663927.3 
#> value: 592567.2 mgc: 151003.6 
#> value: 360201.1 mgc: 150414 
#> value: 154203.4 mgc: 141068 
#> value: 76910.59 mgc: 73580.35 
#> value: 49630.06 mgc: 7802.055 
#> value: 24883.5 mgc: 7316.349 
#> value: 22031.42 mgc: 7907.618 
#> value: 9340.313 mgc: 2979.117 
#> value: 5899.55 mgc: 1764.597 
#> value: 3992.803 mgc: 502.0657 
#> value: 2777.583 mgc: 660.3848 
#> value: 2308.115 mgc: 322.1401 
#> value: 1688.901 mgc: 138.9851 
#> value: 1525.79 mgc: 53.88632 
#> value: 1471.817 mgc: 7.037167 
#> value: 1445.501 mgc: 44.51218 
#> value: 1425.355 mgc: 716.6695 
#> value: 1425.291 mgc: 18.56371 
#> value: 1409.477 mgc: 28.93731 
#> value: 1396.167 mgc: 53.65045 
#> value: 1384.905 mgc: 229.0825 
#> value: 1384.831 mgc: 7.669506 
#> value: 1375.103 mgc: 79.23473 
#> value: 1367.607 mgc: 3.710414 
#> value: 1358.486 mgc: 98.72886 
#> value: 1350.68 mgc: 64.2792 
#> value: 1343.844 mgc: 9.320342 
#> value: 1337.793 mgc: 8.665831 
#> value: 1332.897 mgc: 230.3057 
#> value: 1332.825 mgc: 3.615517 
#> value: 1328.162 mgc: 62.29262 
#> value: 1327.772 mgc: 2.776267 
#> value: 1322.744 mgc: 105.3076 
#> value: 1321.63 mgc: 6.371372 
#> value: 1317.912 mgc: 16.60089 
#> value: 1314.588 mgc: 453.7064 
#> value: 1314.569 mgc: 12.83325 
#> value: 1311.52 mgc: 13.30162 
#> value: 1308.639 mgc: 11.31073 
#> value: 1306.003 mgc: 11.31955 
#> value: 1303.707 mgc: 7.989277 
#> value: 1301.456 mgc: 12.12161 
#> value: 1299.436 mgc: 36.06842 
#> value: 1297.674 mgc: 51.79927 
#> value: 1297.447 mgc: 1.616824 
#> value: 1295.422 mgc: 7.32246 
#> value: 1293.578 mgc: 14.8865 
#> value: 1291.942 mgc: 11.15654 
#> value: 1290.442 mgc: 12.67298 
#> value: 1289.232 mgc: 5.237035 
#> value: 1288.108 mgc: 1.313293 
#> value: 1287.01 mgc: 25.87518 
#> value: 1285.993 mgc: 2.945176 
#> value: 1285.109 mgc: 13.23904 
#> value: 1284.258 mgc: 8.969918 
#> value: 1283.46 mgc: 144.8279 
#> value: 1283.457 mgc: 2.171368 
#> value: 1282.712 mgc: 5.880579 
#> value: 1282.018 mgc: 18.63309 
#> value: 1281.333 mgc: 16.00913 
#> value: 1280.701 mgc: 3.761866 
#> value: 1280.124 mgc: 7.287855 
#> value: 1279.591 mgc: 18.0435 
#> value: 1279.054 mgc: 6.61677 
#> value: 1278.554 mgc: 15.57544 
#> value: 1278.091 mgc: 9.883806 
#> value: 1277.652 mgc: 27.47137 
#> value: 1277.648 mgc: 1.387491 
#> value: 1277.237 mgc: 5.397217 
#> value: 1276.847 mgc: 24.92788 
#> value: 1276.545 mgc: 1.323897 
#> value: 1271.986 mgc: 33.73775 
#> value: 1270.26 mgc: 1.999734 
#> value: 1270.144 mgc: 13.15872 
#> value: 1270.14 mgc: 0.7204035 
#> value: 1270.045 mgc: 3.000105 
#> value: 1269.955 mgc: 13.64413 
#> value: 1269.955 mgc: 0.5588656 
#> value: 1269.879 mgc: 12.71953 
#> value: 1269.859 mgc: 0.4830983 
#> value: 1269.744 mgc: 80.53164 
#> value: 1269.743 mgc: 1.947211 
#> value: 1269.636 mgc: 0.4455192 
#> value: 1269.532 mgc: 48.30398 
#> value: 1269.53 mgc: 1.530554 
#> value: 1269.435 mgc: 4.408405 
#> value: 1269.336 mgc: 3.524596 
#> value: 1269.244 mgc: 3.964306 
#> value: 1269.138 mgc: 21.49671 
#> value: 1269.128 mgc: 1.034898 
#> value: 1269.04 mgc: 4.932344 
#> value: 1269.038 mgc: 0.1879522 
#> value: 1268.781 mgc: 13.78388 
#> value: 1268.607 mgc: 0.8582689 
#> value: 1268.574 mgc: 3.428539 
#> value: 1268.547 mgc: 10.65905 
#> value: 1268.547 mgc: 0.4780089 
#> value: 1268.521 mgc: 2.412192 
#> value: 1268.497 mgc: 0.3939968 
#> value: 1268.473 mgc: 8.043495 
#> value: 1268.47 mgc: 0.2732363 
#> value: 1268.44 mgc: 18.37112 
#> value: 1268.434 mgc: 1.142214 
#> value: 1268.412 mgc: 2.029332 
#> value: 1268.388 mgc: 3.180423 
#> value: 1268.365 mgc: 1.242234 
#> value: 1268.345 mgc: 0.5126696 
#> value: 1268.324 mgc: 0.4391103 
#> value: 1268.305 mgc: 4.543721 
#> value: 1268.292 mgc: 0.230684 
#> value: 1268.228 mgc: 50.42631 
#> value: 1268.228 mgc: 2.551755 
#> value: 1268.172 mgc: 0.3338634 
#> value: 1268.126 mgc: 6.026182 
#> value: 1268.088 mgc: 0.2713439 
#> value: 1268.035 mgc: 10.39659 
#> value: 1268.025 mgc: 0.4893303 
#> value: 1267.988 mgc: 1.080293 
#> value: 1267.959 mgc: 1.919685 
#> value: 1267.932 mgc: 2.455357 
#> value: 1267.911 mgc: 10.5314 
#> value: 1267.903 mgc: 0.5010142 
#> value: 1267.89 mgc: 11.15183 
#> value: 1267.889 mgc: 0.4864803 
#> value: 1267.878 mgc: 1.825475 
#> value: 1267.867 mgc: 17.76147 
#> value: 1267.866 mgc: 1.645518 
#> value: 1267.857 mgc: 1.885765 
#> value: 1267.848 mgc: 0.1500001 
#> value: 1267.836 mgc: 0.9098671 
#> value: 1267.824 mgc: 0.190068 
#> value: 1267.813 mgc: 9.423086 
#> value: 1267.812 mgc: 0.458711 
#> value: 1267.802 mgc: 0.9599396 
#> value: 1267.792 mgc: 22.61884 
#> value: 1267.792 mgc: 0.2377696 
#> value: 1267.783 mgc: 0.6577689 
#> value: 1267.774 mgc: 0.33732 
#> value: 1267.765 mgc: 4.81276 
#> value: 1267.765 mgc: 0.05260067 
#> value: 1267.754 mgc: 0.4577919 
#> value: 1267.745 mgc: 2.220954 
#> value: 1267.745 mgc: 0.0641027 
#> value: 1267.734 mgc: 1.946944 
#> value: 1267.724 mgc: 0.4196361 
#> value: 1267.716 mgc: 6.087933 
#> value: 1267.716 mgc: 0.6407287 
#> value: 1267.709 mgc: 13.84851 
#> value: 1267.709 mgc: 0.7590716 
#> value: 1267.701 mgc: 3.325148 
#> value: 1267.694 mgc: 0.1835423 
#> value: 1267.684 mgc: 0.7930306 
#> value: 1267.674 mgc: 0.4343502 
#> value: 1267.664 mgc: 33.18544 
#> value: 1267.664 mgc: 0.5055578 
#> value: 1267.655 mgc: 1.541604 
#> value: 1267.647 mgc: 0.5222939 
#> value: 1267.64 mgc: 0.3163822 
#> value: 1267.634 mgc: 1.425528 
#> value: 1267.628 mgc: 0.1411388 
#> value: 1267.622 mgc: 0.4182761 
#> value: 1267.617 mgc: 0.5258442 
#> value: 1267.613 mgc: 2.854774 
#> value: 1267.609 mgc: 0.1980784 
#> value: 1267.597 mgc: 0.8861654 
#> value: 1267.587 mgc: 2.807737 
#> value: 1267.586 mgc: 0.1681381 
#> value: 1267.583 mgc: 3.569363 
#> value: 1267.583 mgc: 0.1180857 
#> value: 1267.578 mgc: 0.321002 
#> value: 1267.573 mgc: 0.5474265 
#> value: 1267.569 mgc: 43.35063 
#> value: 1267.569 mgc: 0.1704768 
#> value: 1267.566 mgc: 0.9920958 
#> value: 1267.563 mgc: 3.101464 
#> value: 1267.563 mgc: 0.2371388 
#> value: 1267.561 mgc: 0.5493272 
#> value: 1267.558 mgc: 0.6097116 
#> value: 1267.556 mgc: 0.4706305 
#> value: 1267.554 mgc: 9.366882 
#> value: 1267.554 mgc: 0.5738787 
#> value: 1267.552 mgc: 0.1528771 
#> value: 1267.55 mgc: 0.1335229 
#> value: 1267.548 mgc: 0.5279556 
#> value: 1267.546 mgc: 6.489676 
#> value: 1267.546 mgc: 0.5906358 
#> value: 1267.544 mgc: 6.970496 
#> value: 1267.544 mgc: 0.516262 
#> value: 1267.542 mgc: 0.2784254 
#> value: 1267.539 mgc: 1.161649 
#> value: 1267.537 mgc: 3.88863 
#> value: 1267.537 mgc: 0.09139713 
#> value: 1267.535 mgc: 0.9254993 
#> value: 1267.533 mgc: 0.1457956 
#> value: 1267.531 mgc: 0.4526357 
#> value: 1267.529 mgc: 0.2739218 
#> value: 1267.526 mgc: 0.1988614 
#> value: 1267.524 mgc: 1.26461 
#> value: 1267.523 mgc: 0.6468486 
#> value: 1267.521 mgc: 0.3147832 
#> value: 1267.519 mgc: 1.359603 
#> value: 1267.517 mgc: 0.06082461 
#> value: 1267.51 mgc: 3.944497 
#> value: 1267.506 mgc: 0.1757769 
#> value: 1267.503 mgc: 8.378352 
#> value: 1267.503 mgc: 0.7946084 
#> value: 1267.5 mgc: 0.2005704 
#> value: 1267.498 mgc: 2.766794 
#> value: 1267.498 mgc: 0.06094397 
#> value: 1267.495 mgc: 2.690274 
#> value: 1267.495 mgc: 0.1045404 
#> value: 1267.493 mgc: 0.2379606 
#> value: 1267.491 mgc: 0.9025408 
#> value: 1267.488 mgc: 0.7355307 
#> value: 1267.486 mgc: 0.1450908 
#> value: 1267.484 mgc: 3.370561 
#> value: 1267.484 mgc: 0.2565902 
#> value: 1267.482 mgc: 1.017969 
#> value: 1267.479 mgc: 0.2123067 
#> value: 1267.477 mgc: 0.7938712 
#> value: 1267.475 mgc: 0.3612189 
#> value: 1267.473 mgc: 1.293369 
#> value: 1267.471 mgc: 1.875097 
#> value: 1267.471 mgc: 0.04682026 
#> value: 1267.469 mgc: 0.5072438 
#> value: 1267.467 mgc: 0.07069927 
#> value: 1267.465 mgc: 1.088561 
#> value: 1267.463 mgc: 0.04207041 
#> value: 1267.443 mgc: 36.0518 
#> value: 1267.443 mgc: 2.730045 
#> value: 1267.427 mgc: 0.08506821 
#> value: 1267.425 mgc: 1.3402 
#> value: 1267.425 mgc: 0.08655548 
#> value: 1267.422 mgc: 0.7728272 
#> value: 1267.42 mgc: 0.04267999 
#> value: 1267.418 mgc: 0.6294577 
#> value: 1267.416 mgc: 0.04979805 
#> value: 1267.406 mgc: 1.597226 
#> value: 1267.397 mgc: 0.07300517 
#> value: 1267.395 mgc: 0.9721355 
#> value: 1267.394 mgc: 0.04164901 
#> value: 1267.392 mgc: 2.677317 
#> value: 1267.392 mgc: 0.04846821 
#> value: 1267.39 mgc: 2.27303 
#> value: 1267.39 mgc: 0.1135402 
#> value: 1267.389 mgc: 0.9086769 
#> value: 1267.388 mgc: 0.04185761 
#> value: 1267.386 mgc: 2.716847 
#> value: 1267.386 mgc: 0.1276475 
#> value: 1267.385 mgc: 2.466907 
#> value: 1267.385 mgc: 0.04846515 
#> value: 1267.384 mgc: 1.274862 
#> value: 1267.382 mgc: 0.08011217 
#> value: 1267.381 mgc: 1.121125 
#> value: 1267.38 mgc: 0.07026672 
#> value: 1267.379 mgc: 0.1262073 
#> value: 1267.378 mgc: 0.0659245 
#> value: 1267.378 mgc: 1.379214 
#> value: 1267.378 mgc: 0.06503061 
#> value: 1267.377 mgc: 0.4380098 
#> value: 1267.376 mgc: 0.2328043 
#> value: 1267.375 mgc: 0.2463847 
#> value: 1267.375 mgc: 1.95885 
#> value: 1267.375 mgc: 0.05691374 
#> value: 1267.374 mgc: 0.1390638 
#> value: 1267.373 mgc: 0.5053454 
#> value: 1267.373 mgc: 1.139701 
#> value: 1267.373 mgc: 0.06574726 
#> value: 1267.368 mgc: 39.28754 
#> value: 1267.368 mgc: 1.310191 
#> value: 1267.365 mgc: 0.05627178 
#> value: 1267.364 mgc: 2.406226 
#> value: 1267.363 mgc: 0.1281036 
#> value: 1267.363 mgc: 0.7517985 
#> value: 1267.362 mgc: 0.01913347 
#> value: 1267.361 mgc: 0.2758798 
#> value: 1267.36 mgc: 0.110571 
#> value: 1267.359 mgc: 1.859292 
#> value: 1267.359 mgc: 0.0668008 
#> value: 1267.358 mgc: 0.1274225 
#> value: 1267.357 mgc: 1.179774 
#> value: 1267.357 mgc: 0.05772288 
#> value: 1267.357 mgc: 7.103074 
#> value: 1267.357 mgc: 0.2326094 
#> value: 1267.356 mgc: 0.5548701 
#> value: 1267.355 mgc: 0.169072 
#> value: 1267.355 mgc: 1.384649 
#> value: 1267.355 mgc: 0.04965044 
#> value: 1267.354 mgc: 5.87694 
#> value: 1267.354 mgc: 0.1123729 
#> value: 1267.354 mgc: 1.695678 
#> value: 1267.354 mgc: 0.0464295 
#> value: 1267.353 mgc: 1.505321 
#> value: 1267.353 mgc: 0.04147998 
#> value: 1267.352 mgc: 0.8458953 
#> value: 1267.352 mgc: 0.02957149 
#> value: 1267.351 mgc: 5.420434 
#> value: 1267.351 mgc: 0.3100229 
#> value: 1267.351 mgc: 0.4006729 
#> value: 1267.35 mgc: 0.4740931 
#> value: 1267.349 mgc: 0.05979546 
#> value: 1267.349 mgc: 0.5049684 
#> value: 1267.348 mgc: 0.7675929 
#> value: 1267.347 mgc: 3.224137 
#> value: 1267.347 mgc: 0.2078309 
#> value: 1267.347 mgc: 0.1154254 
#> value: 1267.346 mgc: 0.1853998 
#> value: 1267.346 mgc: 0.5521956 
#> value: 1267.345 mgc: 0.8802601 
#> value: 1267.344 mgc: 0.05039581 
#> value: 1267.344 mgc: 0.5070966 
#> value: 1267.343 mgc: 0.03091531 
#> value: 1267.343 mgc: 0.2907092 
#> value: 1267.342 mgc: 0.7324482 
#> value: 1267.341 mgc: 36.00393 
#> value: 1267.341 mgc: 2.576899 
#> value: 1267.341 mgc: 0.1818308 
#> value: 1267.34 mgc: 0.03196487 
#> value: 1267.34 mgc: 0.4306464 
#> value: 1267.339 mgc: 1.678724 
#> value: 1267.339 mgc: 0.07937939 
#> value: 1267.339 mgc: 1.965848 
#> value: 1267.339 mgc: 0.07790287 
#> value: 1267.338 mgc: 0.1415751 
#> value: 1267.338 mgc: 15.18634 
#> value: 1267.338 mgc: 0.2223445 
#> value: 1267.337 mgc: 0.178543 
#> value: 1267.337 mgc: 3.361097 
#> value: 1267.337 mgc: 0.0543379 
#> value: 1267.336 mgc: 0.3863202 
#> value: 1267.336 mgc: 0.2066731 
#> value: 1267.335 mgc: 2.723113 
#> value: 1267.335 mgc: 0.06216811 
#> value: 1267.335 mgc: 1.154198 
#> value: 1267.335 mgc: 0.05110498 
#> value: 1267.334 mgc: 0.1305264 
#> value: 1267.334 mgc: 0.4646031 
#> value: 1267.333 mgc: 0.1013669 
#> value: 1267.333 mgc: 0.3248838 
#> value: 1267.332 mgc: 1.332741 
#> value: 1267.332 mgc: 0.0792638 
#> value: 1267.332 mgc: 0.2469052 
#> value: 1267.332 mgc: 0.08218199 
#> value: 1267.331 mgc: 0.2763074 
#> value: 1267.331 mgc: 2.456232 
#> value: 1267.331 mgc: 0.1203033 
#> value: 1267.33 mgc: 0.4327469 
#> value: 1267.33 mgc: 0.2734012 
#> value: 1267.33 mgc: 2.065591 
#> value: 1267.33 mgc: 0.1623543 
#> value: 1267.329 mgc: 0.2492188 
#> value: 1267.329 mgc: 1.022884 
#> value: 1267.329 mgc: 0.06033905 
#> value: 1267.328 mgc: 4.735754 
#> value: 1267.328 mgc: 0.1191048 
#> value: 1267.328 mgc: 0.2145898 
#> value: 1267.327 mgc: 3.315181 
#> value: 1267.327 mgc: 0.1034758 
#> value: 1267.327 mgc: 0.09816467 
#> value: 1267.326 mgc: 0.1833721 
#> value: 1267.326 mgc: 0.09216511 
#> value: 1267.326 mgc: 0.05990282 
#> value: 1267.325 mgc: 0.7880752 
#> value: 1267.325 mgc: 0.03315138 
#> value: 1267.325 mgc: 0.2601117 
#> value: 1267.324 mgc: 0.3478648 
#> value: 1267.324 mgc: 5.224672 
#> value: 1267.324 mgc: 0.3186263 
#> value: 1267.324 mgc: 0.05951982 
#> value: 1267.323 mgc: 0.1619244 
#> value: 1267.323 mgc: 0.1021333 
#> value: 1267.322 mgc: 1.75757 
#> value: 1267.322 mgc: 0.0589379 
#> value: 1267.321 mgc: 0.1543897 
#> value: 1267.321 mgc: 1.595594 
#> value: 1267.321 mgc: 0.04362338 
#> value: 1267.321 mgc: 0.09973457 
#> value: 1267.32 mgc: 3.320548 
#> value: 1267.32 mgc: 0.1465666 
#> value: 1267.32 mgc: 0.275159 
#> value: 1267.32 mgc: 2.444444 
#> value: 1267.319 mgc: 0.169374 
#> value: 1267.319 mgc: 0.8856646 
#> value: 1267.319 mgc: 0.03836001 
#> value: 1267.319 mgc: 1.415484 
#> value: 1267.318 mgc: 0.05771761 
#> value: 1267.318 mgc: 0.1841386 
#> value: 1267.318 mgc: 3.577176 
#> value: 1267.318 mgc: 0.06269276 
#> value: 1267.317 mgc: 0.8460506 
#> value: 1267.317 mgc: 0.05341851 
#> value: 1267.317 mgc: 1.017335 
#> value: 1267.317 mgc: 0.03925536 
#> value: 1267.316 mgc: 0.1439491 
#> value: 1267.316 mgc: 0.1193836 
#> value: 1267.316 mgc: 0.06187685 
#> value: 1267.315 mgc: 0.5210019 
#> value: 1267.315 mgc: 0.05018722 
#> value: 1267.315 mgc: 0.0961183 
#> value: 1267.314 mgc: 0.205766 
#> value: 1267.314 mgc: 3.55181 
#> value: 1267.314 mgc: 0.1620508 
#> value: 1267.314 mgc: 0.2633271 
#> value: 1267.313 mgc: 0.6460342 
#> value: 1267.313 mgc: 0.04082085 
#> value: 1267.313 mgc: 2.848909 
#> value: 1267.313 mgc: 0.05326432 
#> value: 1267.312 mgc: 0.06965334 
#> value: 1267.312 mgc: 0.07558166 
#> value: 1267.312 mgc: 0.6273978 
#> value: 1267.312 mgc: 0.02927192 
#> value: 1267.312 mgc: 0.1776044 
#> value: 1267.311 mgc: 0.2524507 
#> value: 1267.311 mgc: 0.5714615 
#> value: 1267.311 mgc: 0.03634663 
#> value: 1267.31 mgc: 6.569855 
#> value: 1267.31 mgc: 0.1836578 
#> value: 1267.31 mgc: 0.4845608 
#> value: 1267.31 mgc: 0.2200542 
#> value: 1267.309 mgc: 0.08006179 
#> value: 1267.309 mgc: 2.06136 
#> value: 1267.309 mgc: 0.06609773 
#> value: 1267.309 mgc: 0.1448974 
#> value: 1267.309 mgc: 0.06329128 
#> value: 1267.308 mgc: 1.120575 
#> value: 1267.308 mgc: 0.055573 
#> value: 1267.308 mgc: 0.5240807 
#> value: 1267.308 mgc: 0.8907478 
#> value: 1267.308 mgc: 0.04499826 
#> value: 1267.308 mgc: 0.503043 
#> value: 1267.307 mgc: 0.02490453 
#> value: 1267.307 mgc: 0.09966598 
#> value: 1267.306 mgc: 0.02953167 
#> value: 1267.306 mgc: 0.1516249 
#> value: 1267.306 mgc: 0.03704595 
#> value: 1267.305 mgc: 1.26063 
#> value: 1267.305 mgc: 0.04350208 
#> value: 1267.305 mgc: 0.4473574 
#> value: 1267.304 mgc: 1.945948 
#> value: 1267.304 mgc: 0.06838389 
#> value: 1267.304 mgc: 0.8137347 
#> value: 1267.304 mgc: 0.04319604 
#> value: 1267.304 mgc: 0.3106419 
#> value: 1267.304 mgc: 0.1588633 
#> value: 1267.303 mgc: 0.1902695 
#> value: 1267.303 mgc: 0.2218362 
#> value: 1267.303 mgc: 0.1798331 
#> value: 1267.303 mgc: 0.1387916 
#> value: 1267.302 mgc: 17.63813 
#> value: 1267.302 mgc: 0.6930635 
#> value: 1267.302 mgc: 0.03094613 
#> value: 1267.299 mgc: 6.404396 
#> value: 1267.298 mgc: 0.4322432 
#> value: 1267.296 mgc: 0.03137295 
#> value: 1267.294 mgc: 0.6520718 
#> value: 1267.292 mgc: 0.04491016 
#> value: 1267.29 mgc: 3.008786 
#> value: 1267.29 mgc: 0.1438102 
#> value: 1267.288 mgc: 0.0311361 
#> value: 1267.287 mgc: 0.1132093 
#> value: 1267.286 mgc: 0.1403634 
#> value: 1267.286 mgc: 0.03119828 
#> value: 1267.286 mgc: 0.8872719 
#> value: 1267.286 mgc: 0.04546892 
#> value: 1267.286 mgc: 0.04113355 
#> value: 1267.285 mgc: 0.1791689 
#> value: 1267.285 mgc: 1.114783 
#> value: 1267.285 mgc: 0.04825858 
#> value: 1267.285 mgc: 0.2973368 
#> value: 1267.285 mgc: 0.04309803 
#> value: 1267.285 mgc: 0.1115932 
#> value: 1267.284 mgc: 0.1243357 
#> value: 1267.284 mgc: 4.801192 
#> value: 1267.284 mgc: 0.4416367 
#> value: 1267.284 mgc: 0.6995811 
#> value: 1267.284 mgc: 0.03277868 
#> value: 1267.284 mgc: 6.548815 
#> value: 1267.284 mgc: 0.470584 
#> value: 1267.284 mgc: 0.0219272 
#> Not improving much - will try early exit...PD hess?: TRUE 

# Compare the estimates and speed
matplot( cbind(opt1$par, opt2$par), type = "l", col = c("black","blue","red"), lty = "solid" )

c(opt1$runtime, opt2$runtime)
#> Time differences in secs
#> [1]  5.26121 39.54170
```
