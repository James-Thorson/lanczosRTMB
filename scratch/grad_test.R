
library(lanczosRTMB)
library(Matrix)

##############
# Run script for testing locally
##############

# Simulate lognormal-gamma process
set.seed(123)
library(RTMB)
n = 30
n_sum = 3
u = 0 + rnorm(n)
y = rgamma( n, shape = 1/0.5^2, scale = exp(u) * 0.5^2 )

# Fit as GLMM
what = "jnll"
nll = function(p){
  sumexpu = sum(exp(p$u[seq_len(n_sum)]))
  nll1 = dnorm(p$u, mean=p$mu, sd=exp(p$logsd), log=TRUE)
  nll2 = dgamma(y, shape = 1/exp(2*p$logcv), scale = exp(p$u) * exp(2*p$logcv), log=TRUE)
  jnll = -1 * ( sum(nll1) + sum(nll2) )
  return(jnll)
}

# Build
obj = lanczos_MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u", k = 10 )
opt = nlminb( obj$par, obj$fn )

# Compare with RTMB
obj2 = MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u", silent = TRUE )
opt2 = nlminb( obj2$par, obj2$fn, obj2$gr )
opt$par - opt2$par

######## Gradient of log-det
# Re-do as penalized likelihood
#newmap = list(mu = factor(NA), logsd = factor(NA), logcv = factor(NA))
#pen = RTMB::MakeADFun( nll, parameters, map = newmap )
#Hq = make_Hq(pen)
#Hq = obj$Hq
#logdet = lanczos_logdet( Hq, k = 10, m = 3, n = length(pen$par), return_extra = FALSE )

#attr(Hq,"env")


###############
# Load objects for lanczos_MakeADFun
###############

func = nll
parameters = obj2$env$parList()
random = "u"
k = 10
profile = NULL
m = 3
seed = 123
do_grad = FALSE

# contents of lanczos_MakeADFun

  #
  env <- new.env(parent = emptyenv())
  fixed = setdiff( names(parameters), c(random,profile) )
  cmb <- function(f, ...) function(p) f(p, ...) ## Helper to make closure

  # DataEval updates
  env$x = unlist(parameters)
  names(env$x) = unlist(sapply( seq_along(parameters), \(i) rep(names(parameters)[i], length(parameters[[i]])) ))
  fetch_x = function() env$x
  #fetch_fixed_vec = function() env$fixed_vec

  # jnll w.r.t. random effects (`parnames = random`) or
  # random and profile (`parnames = c(random,profile)`), conditional on fixed effects
  jnll_vec = function( vec, func, parnames ){
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")
    x = DataEval(fetch_x)
    x[which(names(x) %in% parnames)] = vec
    parlist = parameters
    for(i in seq_along(parlist)){
      parlist[[i]][] = x[which(names(x)==names(parameters)[i])]
    }
    func(parlist )
  }


  # Get tape w.r.t. profile and random effects for optimizing inner problem
  tape_pu = MakeTape(
    f = cmb( jnll_vec, func = func, parnames = c(random, profile) ),
    x = unlist(parameters[names(parameters) %in% c(random, profile)])
  )
  tape_pu$simplify()
  tape_pu$reorder()

  if( isTRUE(TRUE) ){
    # get tape w.r.t. fixed given random effects
    tape_v = MakeTape(
      f = cmb( jnll_vec, func = func, parnames = fixed ),
      x = unlist(parameters[names(parameters) %in% fixed])
    )
    tape_v$simplify()
    tape_v$reorder()
    grad_v = tape_v$jacfun()

    #Hq_q = make_Hv_phi(
    #  obj_phi = tape_phi,
    #  par_phi = unlist(parameters[names(parameters) %in% fixed]),
    #  u_hat = unlist(parameters[names(parameters) %in% c(random, profile)])
    #)
  }

  # Get gradient of tape w.r.t. fixed and random effects for optimizing inner problem
  dfdpu = tape_pu$jacfun()
  dfdpu$simplify()
  dfdpu$reorder()

  # Make function for Hessian w.r.t. random effects (not profiled vars)
  # Hessian-vector product, for Lanczos log-det of Laplace w.r.t. random effects
  tape_u = MakeTape(
    f = cmb( jnll_vec, func = func, parnames = random ),
    x = unlist(parameters[names(parameters) %in% random])
  )
  tape_u$simplify()
  tape_u$reorder()
  #Hq = make_Hq(
  #  tape = tape_u,
  #  x = unlist(parameters[names(parameters) %in% random])
  #)
  #Hq( unlist(parameters[names(parameters) %in% random]) )

  tape_x = MakeTape( nll, parameters )
  Hq = make_Hq(
    tape = tape_x,
    x = unlist(parameters),
    which_random = which(names(env$x) %in% random)
  )
  Hq( unlist(parameters[names(parameters) %in% random]) )

###########
# new code
###########

#
#tape_x = MakeTape( nll, parameters )
x = env$x
fetch_x = function() env$x
which_random = which( names(env$x) %in% random )
which_fixed = which( names(env$x) %in% fixed )

#Hq = obj$Hq
# attr(Hq,"env")$uhat   _should match_  env$x
L = lanczos_logdet( Hq, k = 10, m = 1, n = length(which_random), return_extra = TRUE )

func = nll
random = "u"
fixed = c("mu", "logsd", "logcv")
u_hat = env$x[which_random]
Q_list = lapply(L$L, \(x) x$Q )
whichpos_list = L$which_pos
n_u = length(which_random)

make_d2fdx2_q <-
function( tape,
          x,
          which_random ){

  # Make environment for passing v without retaping
  #local_env <- new.env(parent = emptyenv())
  #local_env$xhat = x
  #local_env$qprime = rep(0, length(x))
  #fetch_qprime = function() local_env$qprime
  qprime = rep(0, length(x))
  q = qprime[which_random]

  # grad
  dfdx = tape$jacfun()
  #dfdx$simplify()

  # grad * v
  dfdx_qprime = function(xq){
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")
    #qprime = DataEval(fetch_qprime)
    x = xq[seq_along(x)]
    qprime[which_random] = xq[-seq_along(x)]
    dfdx(x) %*% qprime
  }
  tape_dfdx_qprime = MakeTape(
    f = dfdx_qprime,
    x = c(x,q)
  )

  # grad( grad * v )
  d2fdx2_qprime = tape_dfdx_qprime$jacfun()

  # Function to supply v for grad( grad * v )
  #d2fdx2_q <- function(xq) {
  #  x = xq[seq_along(x)]
  #  q = xq[-seq_along(x)]
  #  local_env$qprime[which_random] = q
  #  d2fdx2_qprime$force.update()
  #  return(d2fdx2_qprime(x)[which_random])
  #}
  #d2fdx2_q( rnorm(length(which_random)) )

  ## Users must pass q via env
  #d2fdx2_q <- function(x) {
  #  #x = xq[seq_along(x)]
  #  #q = xq[-seq_along(x)]
  #  #local_env$qprime[which_random] = q
  #  d2fdx2_qprime$force.update()
  #  return(d2fdx2_qprime(x)[which_random])
  #}

  # bundle with environment
  out = d2fdx2_qprime
  #attr(out,"env") = local_env
  return(out)
}

# Build and test
q = rnorm(length(which_random))
#xq = c( env$x, q )
Hq( q )

# Compare with new
d2fdx2_q = make_d2fdx2_q( tape_x, x = x, which_random )
#attr(d2fdx2_q,"env")$qprime[which_random] = q
d2fdx2_q( c(env$x,q) )[seq_along(x)][which_random]

# Try gradient of log-det
logdet_wrt_v <- function(v) {
  "[<-" <- ADoverload("[<-")
  "c"   <- ADoverload("c")

  # Pass v into d2fdx2_q
  x = DataEval(fetch_x)
  x[which_fixed] = v
  #attr(d2fdx2_q,"env")$xhat[which_fixed] = v

  total <- 0
  for(i in seq_along(Q_list)) {
    Q   <- Q_list[[i]]        # plain numeric, stop_grad
    k_i = k

    alpha_k = rep(0, k)
    beta_z = rep(0, k - 1)

    for (j in seq_len(k_i)) {
      q = Q[, j]
      #xq = c(x, q)

      # Pass q via environment
      #attr(d2fdx2_q,"env")$qprime[which_random] = q

      # Pass x explicitly to keep on tape
      #w = d2fdx2_q(x)
      w = d2fdx2_q( c(x,q) )[seq_along(x)][which_random]
      alpha_k[j] = sum(q * w)

      if (j < k_i ) {
        # recompute beta from the residual using saved q vectors
        r = w - alpha_k[j] * q
        beta_z[j] <- sqrt(sum(r * r))
      }
   }

    Tri = diag(alpha_k)
    Tri[cbind(2:k,1:(k-1))] = beta_z
    Tri[cbind(1:(k-1),2:k)] = beta_z
    eig  <- eigen(Tri, symmetric = TRUE)
    total <- total + sum(log(eig$values[whichpos_list[[i]]]) * eig$vectors[1,whichpos_list[[i]]]^2) * n_u
  }

  total / m
}
logdet_wrt_v( env$x[which_fixed] )

# Now tape the whole thing w.r.t. phi
tape_logdet <- MakeTape(logdet_wrt_v, x = env$x[which_fixed] )
#tape_logdet$simplify()
#tape_logdet$reorder()
tape_logdet(unlist(parameters[fixed]))

# Gradient function
grad_logdet <- tape_logdet$jacfun()
#grad_logdet$simplify()
#grad_logdet$reorder()

#########
# Compare log-det
#########

# TMB log-det at MLE
determinant(obj2$env$spHess(par = obj2$env$last.par.best, random=TRUE))$modulus

# Lanczos log-det at MLE
Hq = make_Hq(
  tape = tape_x,
  x = obj2$env$last.par.best,
  which_random = which(names(env$x) %in% random)
)
lanczos_logdet( Hq, k = 10, m = 1, n = sum(names(env$x) %in% random) )

#
logdet_wrt_v( unlist(parameters[fixed]) )

#get_logdet = function(v, ... ){
#  attr(Hq,"env")[which_fixed]
#  out = lanczos_logdet(Hq, ...)
#}
#get_logdet( v = )

#########
# Compare grad and FD Test ... seems right
#########

library(numDeriv)
v = obj$par + 1
grad_logdet( v )
grad( logdet_wrt_v, v )

#

#########
# Compare grad and FD Test for total ... not right
#########

v = obj$par
grad_logdet( v ) + grad_v( v )
# Gradient using FD with exact or Lanczos
grad( obj$fn, v )
grad( obj2$fn, v )


#
## Test again
#attr(d2fdx2_q,"env")$xhat[]
#
##
#which_fixed = which( names(env$x) %in% fixed)
#v = opt$par
#Hq <- function(q) {
#  env$q = q
#  d2fdu2_q$force.update()
#  return(d2fdu2_q(env$uhat))
#}

# Debug NaN
DF = grad_logdet$data.frame()
first_bad = min( which(is.nan(DF$Value)) )
DF[ first_bad + seq(-100,100), ]


###################
# Experiment
###################

# Re-do as penalized likelihood
newmap = list(mu = factor(NA), logsd = factor(NA), logcv = factor(NA))
pen = RTMB::MakeADFun( nll, parameters, map = newmap )
Hq = make_Hq(pen)
L = lanczos_logdet( Hq, k = 10, m = 3, n = length(pen$par), return_extra = TRUE )

func = nll
#parameters = parameters
random = "u"
fixed = c("mu", "logsd", "logcv")
u_hat = env$x[which_random]
Q_list = lapply(L$L, \(x) x$Q )
#L_list = lapply(L$L, \(x) x[c("alpha","beta")] )
n_u = length(parameters$u)

  logdet_wrt_v <- function(v) {
    "[<-" <- ADoverload("[<-")
    "c"   <- ADoverload("c")
    u = DataEval(fetch_u)          # stop_grad on u
    x[which_random] = u
    x[which_fixed] = v

    total <- 0
    for (i in seq_len(m)) {
      Q   <- env$Q_list[[i]]        # plain numeric, stop_grad
      #L   <- env$L_list[[i]]
      #k_i <- L$m
      k_i = k

      alphas <- numeric(k_i)
      betas  <- numeric(k_i - 1)

      for (j in seq_len(k_i)) {
        q_j      <- Q[, j]                    # plain numeric
        grad_all <- dfdx(x)              # live w.r.t. phi
        grad_u   <- grad_all[which_random]    # u-block
        alphas[j] <- sum(q_j * grad_u)        # q' H_uu q via grad_u at x_eval

        if (j > 1) {
          # recompute beta from the residual using saved q vectors
          q_prev <- Q[, j-1]
          r <- grad_u - alphas[j] * q_j - betas[j-1] * q_prev
          betas[j-1] <- sqrt(sum(r * r))
        }
      }

      #Tri  <- tridiag(alphas, betas)
      Tri = diag(alphas)
      Tri[cbind(2:k,1:(k-1))] = betas
      Tri[cbind(1:(k-1),2:k)] = betas
      eig  <- eigen(Tri, symmetric = TRUE)
      #wpos <- which(eig$values > 0)
      total <- total +
        sum(log(eig$values) * eig$vectors[1, ]^2) * n_u
    }

    total / m
  }

  # Now tape the whole thing w.r.t. phi
  tape_logdet <- MakeTape(logdet_wrt_v, x = unlist(parameters[fixed]))
  tape_logdet$simplify()
  tape_logdet$reorder()

  # Gradient function
  grad_logdet <- tape_logdet$jacfun()
  grad_logdet$simplify()
  grad_logdet$reorder()


