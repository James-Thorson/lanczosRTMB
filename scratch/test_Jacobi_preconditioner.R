
# pak::pak("james-thorson/lanczosRTMB")

library(Matrix)
library(lanczosRTMB)

# Settings
set.seed(123)
nx = 20
ny = 100
rho = 0.9
logmu = 5

# Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
Px = bandSparse( n = nx, k = c(-1,1), diagonals = list(rep(0.5,nx),rep(0.5,nx)) )
Py = bandSparse( n = ny, k = c(-1,1), diagonals = list(rep(0.5,ny),rep(0.5,ny)) )
Ix = Diagonal(n=nx)
Iy = Diagonal(n=ny)
Qx = (Ix - rho * t(Px)) %*% (10 * Ix) %*% (Ix - rho * Px)
Qy = (Iy - rho * t(Py)) %*% (10 * Iy) %*% (Iy - rho * Py)
Q = kronecker( Qy, Qx )
x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
y = rpois( nx*ny, lambda = exp(logmu + x) )
which_seen = sample( seq_len(nx*ny), size = nx*ny * 0.1, replace = FALSE)
sumy = sum(y)
#sumy = NA
y[-which_seen] = NA

# Function
nll = function(p){
  Qx = (Ix - plogis(p$invlogis_rho) * t(Px)) %*% (exp(2*p$logtau) * Ix) %*% (Ix - plogis(p$invlogis_rho) * Px)
  Qy = (Iy - plogis(p$invlogis_rho) * t(Py)) %*% (exp(2*p$logtau) * Iy) %*% (Iy - plogis(p$invlogis_rho) * Py)
  Q = kronecker( Qy, Qx )
  loglik1 = dgmrf(p$x, Q = Q, log = TRUE)
  loglik2 = sum(dpois(y, exp(p$logmu + p$x), log=TRUE), na.rm=TRUE)
  loglik3 = sum(dnorm(log(sumy), log(sum(exp(p$mu + p$x))), sd = 0.01, log = TRUE), na.rm = TRUE )
  -1 * ( loglik1 + loglik2 + loglik3 )
}
parlist = list( x=rnorm(nx*ny), logtau = 0, invlogis_rho = 0, logmu = 0 )

Hq = make_Hq(
  MakeTape(nll, parlist),
  x0 = parlist,
  which_random = seq_len(nx*ny)
)

dhat = estimate_diag_hutchinson(
  #x = parlist$x,
  Hq = Hq,
  n_probes = 30
)

obj = lanczos_MakeADFun(
  nll,
  parlist,
  k = 40,
  random = "x",
  make_gr = FALSE
)
env = obj$env

inner_opt = newton_CG(
  par = env$x[ names(env$x) %in% "x" ],
  fn = env$tape_pu,
  gr = env$grad_pu,
  Hq = env$Hq_pu,
  diagnostics = TRUE,
  jacobi = 30
)
sum(inner_opt$diag$CG_iter)

start_time = Sys.time()
opt = nlminb(
  obj$par,
  \(x) obj$fn(x, jacobi = 0),
  #obj$gr,
  control = list(trace = 1)
)
( run_time = Sys.time() - start_time )
