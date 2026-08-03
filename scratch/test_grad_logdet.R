
library(Matrix)
library(RTMB)
library(lanczosRTMB)
library(numDeriv)
library(memprof)
source( R'(C:\Users\jtuth\Documents\GitHub\lanczosRTMB\R\lanczos.R)' )

set.seed(123)
n = 100
n_sum = 3

# log-linked normal deviates
u = 0 + rnorm(n)

# Conditional gamma process
y = rgamma( n, shape = 1/0.5^2, scale = exp(u) * 0.5^2 )

nll = function(p){
  sumexpu = sum(exp(p$u[seq_len(n_sum)]))
  ADREPORT( sumexpu )
  REPORT( sumexpu )
  nll1 = dnorm(p$u, mean=p$mu, sd=exp(p$logsd), log=TRUE)
  nll2 = dgamma(
    y,
    shape = 1/exp(2*p$logcv),
    scale = exp(p$u) * exp(2*p$logcv),
    log=TRUE
  )
  jnll = -1 * ( sum(nll1) + sum(nll2) )
  if(what == "jnll") return(jnll)
  if(what == "sumexpu") return(sumexpu)
  if(what == "biascorr") return(jnll + p$eps * sumexpu)
}

# Starting list of parameters
parlist = list(u=u*0, mu = 0, logsd = 0, logcv = 0, eps = 0)

what = "jnll"
obj = lanczos_MakeADFun(
  nll,
  parlist,
  random = "u",
  map = list(eps = factor(NA)),
  k = 10
)

#Hq_u = obj$env$Hq_u
L_extra = lanczos_logdet(
  Hq = obj$env$Hq_u,
  x = obj$env$x[ names(obj$env$x) %in% "u" ],
  k = 10,
  m = 3,
  return_extra = TRUE,
  seed = 123
)

logdet = lanczos_logdet(
  Hq = obj$env$Hq_u,
  x = obj$env$x[ names(obj$env$x) %in% "u" ],
  k = 10,
  m = 3,
  L_extra = L_extra,
  seed = 123
)

cmb <- function(f, ...) function(p) f(p, ...) ## Helper to make closure
tape_logdet_fixedQ = MakeTape(
  #\(x) mean(cmb(
  #  lanczos_logdet,
  #  Hq = obj$env$Hq_u,
  #  k = 10,
  #  m = 3,
  #  L_extra = L_extra
  #)(x)),
  cmb(
    lanczos_logdet,
    Hq = obj$env$Hq_u,
    k = 10,
    m = 3,
    L_extra = L_extra,
    seed = 123
  ),
  obj$env$x[ names(obj$env$x) %in% "u" ]
)
tape_logdet_fixedQ( obj$env$x[ names(obj$env$x) %in% "u" ] )

# Jacobian
grad_logdet = tape_logdet_fixedQ$jacfun()
grad_logdet( obj$env$x[ names(obj$env$x) %in% "u" ] )

###################
# Compare FD ... not identical even with same seed
###################

jac1 = jacobian(
  cmb(
    lanczos_logdet,
    Hq = obj$env$Hq_u,
    k = 10,
    m = 3,
    L_extra = L_extra,
    seed = 123
  ),
  x = obj$env$x[ names(obj$env$x) %in% "u" ]
)

jac2 = jacobian(
  cmb(
    lanczos_logdet,
    Hq = obj$env$Hq_u,
    k = 10,
    m = 3,
    seed = 123
  ),
  x = obj$env$x[ names(obj$env$x) %in% "u" ]
)
