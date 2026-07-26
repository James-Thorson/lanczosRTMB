
library(Matrix)
library(RTMB)
library(lanczosRTMB)
library(numDeriv)
library(memprof)

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
parlist = list(u=u*0, mu = 0, logsd = 0, logcv = 0)

# Build RTMB object
what = "jnll"
obj = RTMB::MakeADFun(
  nll,
  parlist,
  random = "u",
  profile = "mu",
  silent = TRUE
)

# optimize
opt = nlminb( obj$par, obj$fn, obj$gr )
opt

lan_obj = lanczos_MakeADFun(
  nll,
  parlist,
  random = "u",
  profile = "mu",
  k = 10,
  silent = TRUE,
)

if( FALSE ){
  func = nll
  parameters = parlist
  random = "u"
  k = 10
  profile = "mu"
  m = 3
  method = "newton_CG"
  seed = 123
  make_gr = TRUE
  pu_update = c("implicit", "FD", "exact")[1]
  silent = TRUE

}

# optimize
lan_opt = nlminb(
  lan_obj$par,
  \(x) lan_obj$fn(x, orthogonalize = FALSE),
  \(x) lan_obj$gr(x, method = "simple", orthogonalize = FALSE)
)
lan_opt
