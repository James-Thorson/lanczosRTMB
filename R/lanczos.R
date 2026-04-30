


#' @title
#' Assemble tri-diagonal matrix
#'
#' @description
#' Assemble tri-diagonal matrix from alpha and beta from Lanczos method
#'
#' @param Hq function that calculates the product `H %*% q`
#' @param q initial vector used for defining a Kyrlov subspace
#' @param k dimension for Kyrlov subspace
#' @param orthogonalize Whether to do two-pass Gram-Schmidt re-normalization
#'        (much slower)
#' @param tol numerical tolerance for stopping algorithm given that no more terms are identifiable
#'
#' @examples
#' H = diag(exp(rnorm(5)))
#' q = rep(1,5)
#' Hq = function(q) (H %*% q)[,1]
#'
#' L = lanczos(Hq, q, k = 5, ortho = TRUE)
#' T = tridiag(L$alpha, L$beta)
#'
#' # Should match H if and only if L$m = nrow(H)
#' range(L$Q %*% T %*% t(L$Q) - H)
#'
#' @importFrom RTMB MakeADFun sdreport GetTape MakeTape DataEval ADoverload
#' @importFrom Matrix sparseMatrix Diagonal Matrix t
#' @importFrom stats optim rnorm sd
#'
#' @export
lanczos <-
function( Hq,
          q,
          k,
          orthogonalize = FALSE,
          tol = 1e-12 ) {

  n = length(q)
  Q = matrix(0, n, k)
  alpha = numeric(k)
  beta  = numeric(k-1)

  q = q / sqrt(sum(q*q))
  Q[,1] = q
  q_prev = rep(0, n)
  m = k   # actual number of Lanczos steps performed

  for (j in 1:k) {
    w = Hq(q)
    alpha[j] = sum(q * w)

    if (j > 1) {
      w = w - beta[j-1] * q_prev
    }

    if( isTRUE(orthogonalize) ){
      # Full orthogonalize using two-pass Gram-Schmidt
      for (i in 1:j) {
        proj = sum(Q[,i] * w)
        w = w - proj * Q[,i]
      }
      for (i in 1:j) {
        proj = sum(Q[,i] * w)
        w = w - proj * Q[,i]
      }
    }else{
      # Orthogonalize with previous
      # (accumulates errors but faster)
      w = w - alpha[j] * q
    }

    if (j < k) {
      beta[j] = sqrt(sum(w*w))
      if(beta[j] < tol) {
        m = j     # record the number of completed steps
        break     # terminate early
      }
      q_prev = q
      q = w / beta[j]
      Q[,j+1] = q
    }
  }

  list(
    Q = Q[, seq_len(m), drop = FALSE],
    alpha = alpha[seq_len(m)],
    beta = beta[seq_len(m-1)],
    m = m
  )
}

#' @title
#' Make function to calculate H %*% q
#'
#' @description
#' Given a TMB object, make a function that efficiently calculates
#'   `H %*% q` without constructing H itself, and instead using `grad_u( grad_u(f) %** q)`
#'   given function f(x) that returns the negative log-likelihood given `x = u` with fixed `v`.  This can
#'   then be used e.g. in Lanczos methods when H is too large to construct explicitly
#'
#' @param uhat parameter vector `u` used when evaluating `H`
#' @param tape Alternative to specifying `obj`
#'
#' @examples
#' u = rnorm(100)
#' y = rpois(length(u), exp(u))
#' nll = function(p) -1 * ( sum(dnorm(p$u,log=TRUE)) + sum(dpois(y,exp(p$u),log=TRUE)) )
#' obj = RTMB::MakeADFun( nll, list(u=u), silent = TRUE )
#' Hq = make_Hq( GetTape(obj), obj$par )
#' # Confirm
#' q = rnorm( length(obj$par) )
#' all.equal( Hq(q)[1,], (obj$he()%*%q)[,1] )
#'
#' @export
make_Hq <-
function( tape,
          uhat ){

  # Make environment for passing v without retaping
  env <- new.env(parent = emptyenv())
  env$uhat = uhat
  env$q = 0 * uhat
  fetch_q = function() env$q

  # grad
  dfdu = tape$jacfun()
  dfdu$simplify()

  # grad * v
  dfdu_q = function(u){
    q = DataEval(fetch_q)
    dfdu(u) %*% q
  }
  tape_dfdu_q = MakeTape(
    f = dfdu_q,
    x = uhat
  )
  tape_dfdu_q$simplify()
  tape_dfdu_q$reorder()

  # grad( grad * v )
  d2fdu2_q = tape_dfdu_q$jacfun()
  d2fdu2_q$simplify()
  d2fdu2_q$reorder()

  # Function to supply v for grad( grad * v )
  Hq <- function(q) {
    env$q = q
    d2fdu2_q$force.update()
    return(d2fdu2_q(env$uhat))
  }

  # bundle with environment
  out = Hq
  attr(out,"env") = env
  return(out)
}

#' @title
#' Assemble tri-diagonal symmetric matrix
#'
#' @description
#' Assemble tri-diagonal matrix from alpha and beta from Lanczos method
#'
#' @param alpha vector for diagonal
#' @param beta vector for off-diagonal
#'
#' @export
tridiag <-
function( alpha,
          beta ) {
  k = length(alpha)
  Tri = Matrix(0, k, k)
  diag(Tri) = alpha
  for (i in seq_len(k-1)) {
    Tri[i, i+1] = beta[i]
    Tri[i+1, i] = beta[i]
  }
  return(Tri)
}



#' @title
#' Sample from Lanczos method
#'
#' @description
#' Sample from a penalized likelihood model without directly forming the Hessian
#' matrix
#'
#' @inheritParams lanczos
#' @param nsamp number of samples
#'
#' @examples
#' # Simulate lognormal-gamma process
#' set.seed(123)
#' library(RTMB)
#' n = 30
#' n_sum = 3
#' u = 0 + rnorm(n)
#' y = rgamma( n, shape = 1/0.5^2, scale = exp(u) * 0.5^2 )
#'
#' # Fit as GLMM
#' what = "jnll"
#' nll = function(p){
#'   sumexpu = sum(exp(p$u[seq_len(n_sum)]))
#'   ADREPORT( sumexpu )
#'   REPORT( sumexpu )
#'   nll1 = dnorm(p$u, mean=p$mu, sd=exp(p$logsd), log=TRUE)
#'   nll2 = dgamma(y, shape = 1/exp(2*p$logcv), scale = exp(p$u) * exp(2*p$logcv), log=TRUE)
#'   jnll = -1 * ( sum(nll1) + sum(nll2) )
#'   if(what == "jnll") return(jnll)
#'   if(what == "sumexpu") return(sumexpu)
#' }
#' obj = RTMB::MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u", silent = TRUE )
#' opt = nlminb( obj$par, obj$fn, obj$gr )
#' sdrep = sdreport(obj, bias.correct = TRUE )
#' H = obj$env$spHess(par = obj$env$last.par.best, random = TRUE)
#'
#' # Re-do as penalized likelihood
#' newmap = list(mu = factor(NA), logsd = factor(NA), logcv = factor(NA))
#' pen = RTMB::MakeADFun( nll, obj$env$parList(), map = newmap, silent = TRUE )
#' opt_pen = nlminb( pen$par, pen$fn, pen$gr )
#' Hq = make_Hq( GetTape(pen), opt_pen$par )
#'
#' # Gradient-based Lanczos sampling
#' what = "sumexpu"
#' grad = RTMB::MakeADFun( nll, obj$env$parList(), map = newmap, silent = TRUE )$gr( opt_pen$par )[1,]
#' what = "jnll"
#' sample_x = function(n){ x = rnorm(n); return( x / sqrt(sum(x^2))) }
#' samples = lanczos_sample(
#'   Hq = Hq,
#'   q = grad,
#'   k = 30,
#'   n = 1000,
#'   orthogonalize = TRUE
#' )
#'
#' # Samples from Lanczos for parameters that contribute to derived quantity
#' samples = sweep( samples, MARGIN = 1, FUN = "+", STATS = opt$par )
#' apply( samples, MARGIN = 1, FUN = sd )
#'
#' # Samples from full Hessian should approximately match
#' samp2 = t(mvtnorm::rmvnorm( n = 1000, sigma = as.matrix(solve(H)) ))
#' apply( samp2, MARGIN = 1, FUN = sd )
#'
#' # Compare bias-correction
#' sumexpu_z = apply( samples, MARGIN = 2, FUN = \(x) pen$report(x)$sumexpu )
#' mean(sumexpu_z)
#' summary(sdrep)['sumexpu',]
#'
#' @export
lanczos_sample <-
function( Hq,
          q,
          k,
          nsamp = 1,
          orthogonalize = FALSE) {

  norm_q = sqrt(sum(q^2))
  q = q / norm_q

  #
  L = lanczos(Hq, q, k, orthogonalize = orthogonalize)
  T = tridiag(L$alpha, L$beta)
  if( FALSE ){
    L$Q %*% T %*% t(L$Q) # SHOULD MATCH H
  }

  # eigendecomposition of small matrix
  eig = eigen(T, symmetric = TRUE)
  V = eig$vectors
  D = eig$values

  # truncate
  whichpos = which(D>0)
  D = D[whichpos]
  V = V[,whichpos]

  # sample in Krylov space
  xi = matrix( rnorm(length(D)*nsamp), ncol = nsamp )
  y  = V %*% sweep(xi, MARGIN = 1, FUN = "/", STATS = sqrt(D))

  # lift to full space
  z = L$Q %*% y        # norm_v *
  return(z)
}


#' @title
#' Estimate variance using Lanczos method
#'
#' @description
#' Estimate squared standard error (variance) for a parameter or derived quantity
#'    without forming the Hessian matrix directly
#'
#' @inheritParams lanczos
#' @param q vector to use when calculating variance, either an indicator for a single
#'    parameter, or a gradient evaluated at the MLE for a derived quantity
#' @param k can be a vector
#' @param min_spectral_ratio is the ratio of minimum to maximum ratio,
#'    where values < minimum are truncated to minimum
#'    (min_spectral_ratio=0 disables truncation)
#'
#' @export
lanczos_variance <-
function( Hq,
          q,
          k = c(25,30),
          min_spectral_ratio = 1e-10,
          orthogonalize = FALSE ) {

  norm_q = sqrt(sum(q^2))
  q = q / norm_q
  L = lanczos( Hq, q, max(k), orthogonalize = orthogonalize )
  Tri = tridiag(L$alpha, L$beta)

  fn = function( dim ){
    T = Tri[seq_len(dim),seq_len(dim)]
    # invert small matrix safely
    eig = eigen(T, symmetric = TRUE)
    lambda = eig$values

    # Spectral clipping
    eps = min_spectral_ratio * max(lambda)
    lambda = ifelse(lambda<eps, eps, lambda)

    Tinv = eig$vectors %*% diag(1 / lambda) %*% t(eig$vectors)

    # e1' T^{-1} e1 = (1,1) element
    quad = Tinv[1, 1]

    # final variance estimate
    (norm_q^2) * quad
  }
  sapply( X = k, FUN = fn )
}


#' @title
#' Estimate log-determinant likelihood using Lanczos method
#'
#' @description
#' Estimate log-determinant of inner Hessian matrix with a stochastic approximation
#'
#' @inheritParams lanczos
#' @param m number of probe-vectors to use for approximating average and standard
#'        deviation of log-determinant
#' @param n length of parameters (and necessary probe-vector)
#' @param seed if not NULL, then sets the seed.  This is helfpul given that
#'    the Hutchinson probe vectors are randomly sampled, and comparisons have
#'    lower variance using a fixed seed.
#' @param return_extra whether to return probes and other internal constructions.
#'
#' @details
#' For a model with independent random effects, the variance of stochastic
#' trace estimation should be zero
#'
#' @examples
#' # Simulate lognormal-gamma process
#' set.seed(123)
#' library(RTMB)
#' n = 30
#' n_sum = 3
#' u = 0 + rnorm(n)
#' y = rgamma( n, shape = 1/0.5^2, scale = exp(u) * 0.5^2 )
#'
#' # Fit as GLMM
#' what = "jnll"
#' nll = function(p){
#'   sumexpu = sum(exp(p$u[seq_len(n_sum)]))
#'   ADREPORT( sumexpu )
#'   REPORT( sumexpu )
#'   nll1 = dnorm(p$u, mean=p$mu, sd=exp(p$logsd), log=TRUE)
#'   nll2 = dgamma(y, shape = 1/exp(2*p$logcv), scale = exp(p$u) * exp(2*p$logcv), log=TRUE)
#'   jnll = -1 * ( sum(nll1) + sum(nll2) )
#'   if(what == "jnll") return(jnll)
#'   if(what == "sumexpu") return(sumexpu)
#' }
#' obj = RTMB::MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u", silent = TRUE )
#' opt = nlminb( obj$par, obj$fn, obj$gr )
#' sdrep = sdreport(obj, bias.correct = TRUE )
#' H = obj$env$spHess(par = obj$env$last.par.best, random = TRUE)
#'
#' # Re-do as penalized likelihood
#' newmap = list(mu = factor(NA), logsd = factor(NA), logcv = factor(NA))
#' pen = RTMB::MakeADFun( nll, obj$env$parList(), map = newmap, silent = TRUE )
#' opt_pen = nlminb( pen$par, pen$fn, pen$gr )
#'
#' # Compare determinant
#' Hq = make_Hq( GetTape(pen), opt_pen$par )
#' lanczos_logdet( Hq, k = 10, m = 3, n = length(pen$par) )
#' Matrix::determinant( H )
#'
#' @export
lanczos_logdet <-
function( Hq,
          k,
          m,
          n,
          seed = NULL,
          orthogonalize = TRUE,
          return_extra = FALSE ) {

  if( !is.null(seed) ){
    set.seed(seed)
  }
  q_m = matrix( sample( c(-1,1), size = n*m, replace=TRUE), ncol = m )  # Rademacher vector
  logdet = numeric(m)
  which_pos = Tri = eig = L = vector("list", length = m)

  for (mi in 1:m) {
    q = q_m[,mi]
    L[[mi]] = lanczos( Hq = Hq, q = q, k = max(k), orthogonalize = orthogonalize )
    Tri[[mi]] = tridiag(L[[mi]]$alpha, L[[mi]]$beta)
    eig[[mi]] = eigen(Tri[[mi]], symmetric = TRUE)
    which_pos[[mi]] = which( eig[[mi]]$values > 0 )
    log_quad = sum(log(eig[[mi]]$values[which_pos[[mi]]]) * eig[[mi]]$vectors[1, which_pos[[mi]]]^2)
    logdet[mi] = log_quad * n
  }

  #
  if( isFALSE(return_extra) ){
    return( logdet )
  }else{
    return( list(logdet = logdet, q_m = q_m, L = L, Tri = Tri, eig = eig, which_pos = which_pos) )
  }
}

#' @title
#' Estimate log-marginal likelihood using Lanczos method
#'
#' @description
#' Estimate log-marginal likelihood using Laplace approximation, but replacing
#'    exact calculation of log-determinant with a stochastic approximation
#'
#' @details
#' This function is only intended when integrating across all parameters, e.g.,
#' when supplying a penalized likelihood model with fixed effects mapped off at
#' a prior estimate.  For more control over which parameters to estimate, use
#' [lanczos_MakeADFun]
#'
#' @inheritParams lanczos
#' @inheritParams lanczos_logdet
#' @inheritParams make_Hq
#' @param obj output from `TMB::MakeADFun` when using penalized likelihood
#'
#' @examples
#' # Simulate lognormal-gamma process
#' set.seed(123)
#' library(RTMB)
#' n = 30
#' n_sum = 3
#' u = 0 + rnorm(n)
#' y = rgamma( n, shape = 1/0.5^2, scale = exp(u) * 0.5^2 )
#'
#' # Fit as GLMM
#' what = "jnll"
#' nll = function(p){
#'   sumexpu = sum(exp(p$u[seq_len(n_sum)]))
#'   ADREPORT( sumexpu )
#'   REPORT( sumexpu )
#'   nll1 = dnorm(p$u, mean=p$mu, sd=exp(p$logsd), log=TRUE)
#'   nll2 = dgamma(y, shape = 1/exp(2*p$logcv), scale = exp(p$u) * exp(2*p$logcv), log=TRUE)
#'   jnll = -1 * ( sum(nll1) + sum(nll2) )
#'   if(what == "jnll") return(jnll)
#'   if(what == "sumexpu") return(sumexpu)
#' }
#' obj = RTMB::MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u", silent = TRUE )
#' opt = nlminb( obj$par, obj$fn, obj$gr )
#' sdrep = sdreport(obj, bias.correct = TRUE )
#' H = obj$env$spHess(par = obj$env$last.par.best, random = TRUE)
#'
#' # Re-do as penalized likelihood
#' newmap = list(mu = factor(NA), logsd = factor(NA), logcv = factor(NA))
#' pen = RTMB::MakeADFun( nll, obj$env$parList(), map = newmap, silent = TRUE )
#' opt_pen = nlminb( pen$par, pen$fn, pen$gr )
#'
#' # Compare marginal likelihoods
#' lanczos_nll( pen, k = 10, m = 10 )
#' opt$obj
#'
#' @export
lanczos_nll <-
function( obj,
          uhat = obj$env$last.par.best,
          k,
          m,
          seed = NULL ) {

  #
  Hq = make_Hq( GetTape(obj), uhat )
  logdet = lanczos_logdet( Hq, k, m, n = length(uhat), seed = seed )
  nll = obj$env$f(uhat) - (0.5*length(uhat)*log(2*pi)) + 0.5*mean(logdet)
  sd_nll = 0.5 * sd(logdet)
  return( c(nll = nll, sd_nll = sd_nll) )
}

#' @title
#' Approximate log-marginal likelihood using Lanczos method (EXPERIMENTAL)
#'
#' @description
#' Make a function that returns the Lanczos-Laplace approximation for
#'   a user-supplied joint likelihood and designated random effects
#'
#' @inheritParams lanczos_logdet
#' @inheritParams RTMB::MakeADFun
#' @inheritParams TMB::MakeADFun
#'
#' @return
#' An object (list) of class `tinyVAST`. Elements include:
#' \describe{
#' \item{par}{parameter-vector of fixed effects}
#' \item{fn}{function that returns the Lanczos-Laplace approximation given fixed effects}
#' \item{env}{environment of local variables}
#' \item{Hq}{function that returns the Hessian-vector product (for use in debugging)}
#' }
#'
#' @examples
#' # Simulate lognormal-gamma process
#' set.seed(123)
#' library(RTMB)
#' n = 30
#' n_sum = 3
#' u = 0 + rnorm(n)
#' y = rgamma( n, shape = 1/0.5^2, scale = exp(u) * 0.5^2 )
#'
#' # Fit as GLMM
#' what = "jnll"
#' nll = function(p){
#'   sumexpu = sum(exp(p$u[seq_len(n_sum)]))
#'   nll1 = dnorm(p$u, mean=p$mu, sd=exp(p$logsd), log=TRUE)
#'   nll2 = dgamma(y, shape = 1/exp(2*p$logcv), scale = exp(p$u) * exp(2*p$logcv), log=TRUE)
#'   jnll = -1 * ( sum(nll1) + sum(nll2) )
#'   return(jnll)
#' }
#'
#' # Build
#' obj = lanczos_MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u", k = 10 )
#' opt = nlminb( obj$par, obj$fn )
#'
#' # Compare with RTMB
#' obj2 = MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u", silent = TRUE )
#' opt2 = nlminb( obj2$par, obj2$fn, obj2$gr )
#' opt$par - opt2$par
#'
#' @export
lanczos_MakeADFun <-
function( func,
          parameters,
          random,
          k,
          profile = NULL,
          m = 3,
          #do_grad = FALSE,
          seed = 123 ){

  # vectors
  #  p = profile
  #  u = random
  #  v = fixed
  #  x = (p,u,v)

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

  #if( isTRUE(do_grad) ){
  #  # get tape w.r.t. fixed given random effects
  #  tape_v = MakeTape(
  #    f = cmb( jnll_vec, func = func, parnames = fixed ),
  #    x = unlist(parameters[names(parameters) %in% fixed])
  #  )
  #  tape_v$simplify()
  #  tape_v$reorder()
  #  grad_v = tape_v$jacfun()

  #  Hq_q = make_Hv_phi(
  #    obj_phi = tape_phi,
  #    par_phi = unlist(parameters[names(parameters) %in% fixed]),
  #    u_hat = unlist(parameters[names(parameters) %in% c(random, profile)])
  #  )
  #}

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
  Hq = make_Hq(
    tape = tape_u,
    uhat = unlist(parameters[names(parameters) %in% random])
  )

  # Experiment
  #if( FALSE ){
  #  tape_joint = MakeTape(
  #    f = func,
  #    x = parameters
  #  )
  #}

  # Objective function
  nll = function(v){
    # Define fixed effects and assign to global environment
    env$x[which(names(env$x) %in% fixed)] = v
    tape_pu$force.update()
    dfdpu$force.update()

    if( is.null(env$inner_start) ){
      env$puhat = env$x[names(env$x) %in% c(random,profile)]
      env$best = Inf
    }

    # Run inner and assign xhat to global environment
    out = optim(
      par = env$puhat,
      fn = tape_pu,
      gr = dfdpu,
      method = "L-BFGS-B",
      control = list(trace=0, maxit = 1e3, factr = 1e-2)
    )

    # Define profiled effects, and assign to global environment
    which_profile = which( names(out$par) %in% profile )
    if( length(which_profile) > 0 ){
      env$x[which(names(env$x) %in% profile)] = out$par[which_profile]
      tape_u$force.update()
    }

    # Have to assign into env(Hq)$mle to evaluate at right point
    which_random = which( names(out$par) %in% random )
    if( length(which_random) > 0 ){
      attr(Hq,"env")$uhat = out$par[which_random]
      env$logdet1_m = lanczos_logdet(
        Hq = Hq,
        k = k,
        m = m,
        n = length(which_random),
        return_extra = FALSE,
        seed = seed
      )
    }else{
      env$logdet1_m = rep(0, m)
    }

    # Get jnll
    jnll = tape_pu( out$par )

    # Assemble laplace
    #if(isTRUE(do_grad)){
    #  neglogmarglik = jnll - (0.5*length(out$par)*log(2*pi)) + 0.5*mean(env$logdet1_m$logdet)
    #}else{
      neglogmarglik = jnll - (0.5*length(out$par)*log(2*pi)) + 0.5*mean(env$logdet1_m)
    #}

    #
    #if( isTRUE(do_grad) ){
    #  # Gradient of joint likelihood w.r.t. fixed effects
    #  grad_phi$force.update()
    #  grad1 = grad_phi(fixedvec)
    #}

    # Assign best and return
    if( neglogmarglik < env$best ){
      env$puhat = out$par
      env$best = neglogmarglik
    }
    return( neglogmarglik )
  }

  #
  out = list(
    par = unlist(parameters[fixed]),
    fn = nll,
    env = env,
    Hq = Hq
  )
  return(out)
}
