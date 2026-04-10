


#' @title
#' Assemble tri-diagonal matrix
#'
#' @description
#' Assemble tri-diagonal matrix from alpha and beta from Lanczos method
#'
#' @param Hv function that calculates the product `Hv`
#' @param q1 initial vector used for defining a Kyrlov subspace
#' @param k dimension for Kyrlov subspace
#' @param orthogonalize Whether to do two-pass Gram-Schmidt re-normalization
#'        (much slower)
#' @param tol numerical tolerance for stopping algorithm given that no more terms are identifiable
#'
#' @examples
#' H = diag(exp(rnorm(5)))
#' v = rep(1,5)
#' Hv = function(v) (H %*% v)[,1]
#'
#' L = lanczos(Hv, v, k = 5, ortho = TRUE)
#' T = tridiag(L$alpha, L$beta)
#'
#' # Should match H if and only if L$m = nrow(H)
#' range(L$Q %*% T %*% t(L$Q) - H)
#'
#' @importFrom RTMB MakeADFun sdreport GetTape MakeTape DataEval
#' @importFrom Matrix sparseMatrix Diagonal Matrix t
#'
#' @export
lanczos <-
function( Hv,
          q1,
          k,
          orthogonalize = FALSE,
          tol = 1e-12 ) {

  n = length(q1)
  Q = matrix(0, n, k)
  alpha = numeric(k)
  beta  = numeric(k-1)

  q = q1 / sqrt(sum(q1*q1))
  Q[,1] = q
  q_prev = rep(0, n)
  m = k   # actual number of Lanczos steps performed

  for (j in 1:k) {
    w = Hv(q)
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
#' Make function to calculate H %*% v
#'
#' @description
#' Given a TMB object, make a function that efficiently calculates
#'   `H %*% v` without constructing H itself, and instead using `grad_x( grad_x(f) %** v)`
#'   given function f(x) that returns the negative log-likelihood.  This can
#'   then be used e.g. in Lanczos methods when H is too large to construct explicitly
#'
#' @param obj TMB object (output from `TMB::MakeADFun`)
#' @param par parameter vector `p` used when evaluating `H`
#'
#' @examples
#' u = rnorm(100)
#' y = rpois(length(u), exp(u))
#' nll = function(p) -1 * ( sum(dnorm(p$u,log=TRUE)) + sum(dpois(y,exp(p$u),log=TRUE)) )
#' obj = RTMB::MakeADFun( nll, list(u=u) )
#' Hv = make_Hv( obj )
#' # Confirm
#' v = rnorm( length(obj$par) )
#' all.equal( Hv(v)[1,], (obj$he()%*%v)[,1] )
#'
#' @export
make_Hv <-
function( obj,
          par = obj$env$last.par.best,
          tape ){

  # Make environment for passing v without retaping
  env <- new.env(parent = emptyenv())
  env$mle = par
  env$v = 0 * par
  fetch_v = function() env$v

  # grad
  if(missing(tape)) tape = GetTape(obj)
  dfdx = tape$jacfun()
  dfdx$simplify()

  # grad * v
  dfdx_v = function(x){
    v = DataEval(fetch_v)
    dfdx(x) %*% v
  }
  tape_dfdx_v = MakeTape(
    f = dfdx_v,
    x = par
  )
  tape_dfdx_v$simplify()
  tape_dfdx_v$reorder()

  # grad( grad * v )
  d2fdx2_v = tape_dfdx_v$jacfun()
  d2fdx2_v$simplify()
  d2fdx2_v$reorder()

  # Function to supply v for grad( grad * v )
  Hv <- function(v) {
    env$v = v
    d2fdx2_v$force.update()
    return(d2fdx2_v(env$mle))
  }

  # bundle with environment
  out = Hv
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
#' obj = RTMB::MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u" )
#' opt = nlminb( obj$par, obj$fn, obj$gr )
#' sdrep = sdreport(obj, bias.correct = TRUE )
#' H = obj$env$spHess(par = obj$env$last.par.best, random = TRUE)
#'
#' # Re-do as penalized likelihood
#' newmap = list(mu = factor(NA), logsd = factor(NA), logcv = factor(NA))
#' pen = RTMB::MakeADFun( nll, obj$env$parList(), map = newmap )
#' opt_pen = nlminb( pen$par, pen$fn, pen$gr )
#' Hv = make_Hv( pen )
#'
#' # Gradient-based Lanczos sampling
#' what = "sumexpu"
#' grad = RTMB::MakeADFun( nll, obj$env$parList(), map = newmap )$gr( opt$par )[1,]
#' what = "jnll"
#' sample_x = function(n){ x = rnorm(n); return( x / sqrt(sum(x^2))) }
#' samples = lanczos_sample(
#'   Hv = Hv,
#'   v = grad,
#'   k = 30,
#'   n = 1000,
#'   orthogonalize = TRUE
#' )
#'
#' # Samples from Lanczos for parameters that contribute to derived quantity
#' samples = sweep( samples, MARGIN = 1, FUN = "+", STAT = opt$par )
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
function( Hv,
          v,
          k,
          nsamp = 1,
          orthogonalize = FALSE) {

  norm_v = sqrt(sum(v^2))
  v = v / norm_v

  #
  L = lanczos(Hv, v, k, orthogonalize = orthogonalize)
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
  y  = V %*% sweep(xi, MARGIN = 1, FUN = "/", STAT = sqrt(D))

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
#' @param v vector to use when calculating variance, either an indicator for a single
#'    parameter, or a gradient evaluated at the MLE for a derived quantity
#' @param k can be a vector
#' @param min_spectral_ratio is the ratio of minimum to maximum ratio,
#'    where values < minimum are truncated to minimum
#'    (min_spectral_ratio=0 disables truncation)
#'
#' @export
lanczos_variance <-
function( Hv,
          v,
          k = c(25,30),
          min_spectral_ratio = 1e-10,
          orthogonalize = FALSE ) {

  norm_v = sqrt(sum(v^2))
  v = v / norm_v
  L = lanczos( Hv, v, max(k), orthogonalize = orthogonalize )
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
    (norm_v^2) * quad
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
#'
#' @details
#' For a model with independent random effects, the variance of stochastic
#' trace estimation should be zero
#'
#' @export
lanczos_logdet <-
function( Hv,
          k,
          m,
          n,
          seed = NULL,
          orthogonalize = TRUE ) {

  if( !is.null(seed) ){
    set.seed(seed)
  }
  v_m = matrix( sample( c(-1,1), size = n*m, replace=TRUE), ncol = m )  # Rademacher vector
  logdet = numeric(m)

  for (mindex in 1:m) {
    v = v_m[,mindex]
    L = lanczos( Hv = Hv, q1 = v, k = max(k), orthogonalize = orthogonalize )
    Tri = tridiag(L$alpha, L$beta)
    eig = eigen(Tri, symmetric = TRUE)
    which_pos = which( eig$values > 0 )
    log_quad = sum(log(eig$values[which_pos]) * eig$vectors[1, which_pos]^2)
    logdet[mindex] = log_quad * length(v)
  }

  #
  return( logdet )
}

#' @title
#' Estimate log-marginal likelihood using Lanczos method
#'
#' @description
#' Estimate log-marginal likelihood using Laplace approximation, but replacing
#'    exact calculation of log-determinant with a stochastic approximation
#'
#' @inheritParams lanczos
#' @inheritParams lanczos_logdet
#' @inheritParams make_Hv
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
#' obj = RTMB::MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u" )
#' opt = nlminb( obj$par, obj$fn, obj$gr )
#' sdrep = sdreport(obj, bias.correct = TRUE )
#' H = obj$env$spHess(par = obj$env$last.par.best, random = TRUE)
#'
#' # Re-do as penalized likelihood
#' newmap = list(mu = factor(NA), logsd = factor(NA), logcv = factor(NA))
#' pen = RTMB::MakeADFun( nll, obj$env$parList(), map = newmap )
#' opt_pen = nlminb( pen$par, pen$fn, pen$gr )
#' Hv = make_Hv( pen )
#'
#' # Compare determinant
#' Hv = make_Hv(pen)
#' lanczos_logdet( Hv, k = 10, m = 3, n = length(pen$par) )
#' Matrix::determinant( H )
#'
#' # Compare marginal likelihoods
#' lanczos_nll( pen, k = 10, m = 10 )
#' opt$obj
#'
#' @export
lanczos_nll <-
function( obj,
          k,
          m,
          Hv,
          seed = NULL ) {

  if( missing(Hv) ){
    Hv = make_Hv( obj )
  }

  #
  inner_par = obj$env$last.par.best
  logdet = lanczos_logdet( Hv, k, m, n = length(inner_par), seed = seed )
  nll = obj$env$f(inner_par) - (0.5*length(inner_par)*log(2*pi)) + 0.5*mean(logdet)
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
#'
#' @return
#' An object (list) of class `tinyVAST`. Elements include:
#' \describe{
#' \item{par}{parameter-vector of fixed effects}
#' \item{fn}{function that returns the Lanczos-Laplace approximation given fixed effects}
#' \item{env}{environment of local variables}
#' \item{Hv}{function that returns the Hessian-vector product (for use in debugging)}
#' }
#'
#' @export
lanczos_MakeADFun <-
function( func,
          parameters,
          random,
          k,
          profile = NULL,
          m = 3,
          seed = 123 ){

  #
  env <- new.env(parent = emptyenv())
  fixed = setdiff( names(parameters), c(random,profile) )
  cmb <- function(f, ...) function(p) f(p, ...) ## Helper to make closure

  # DataEval updates
  env$par_vec = unlist(parameters)
  names(env$par_vec) = unlist(sapply( seq_along(parameters), \(x) rep(names(parameters)[x], length(parameters[[x]])) ))
  fetch_par_vec = function() env$par_vec
  #fetch_fixed_vec = function() env$fixed_vec

  # jnll w.r.t. random effects (`parnames = random`) or
  # random and profile (`parnames = c(random,profile)`), conditional on fixed effects
  jnll_x = function( xvec, func, parnames ){
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")
    par_vec = DataEval(fetch_par_vec)
    par_vec[which(names(par_vec) %in% parnames)] = xvec
    parlist = parameters
    for(i in seq_along(parlist)){
      parlist[[i]][] = par_vec[which(names(par_vec)==names(parameters)[i])]
    }
    func(parlist )
  }

  # Get tape w.r.t. fixed and random effects for optimizing inner problem
  tape_u = MakeTape(
    f = cmb( jnll_x, func = func, parnames = c(random, profile) ),
    x = unlist(parameters[names(parameters) %in% c(random, profile)])
  )
  tape_u$simplify()
  tape_u$reorder()

  # Get gradient of tape w.r.t. fixed and random effects for optimizing inner problem
  dfdu = tape_u$jacfun()
  dfdu$simplify()
  dfdu$reorder()

  # Make function for Hessian w.r.t. random effects (not profiled vars)
  # Hessian-vector product, for Lanczos log-det of Laplace w.r.t. random effects
  tape_random = MakeTape(
    f = cmb( jnll_x, func = func, parnames = random ),
    x = unlist(parameters[names(parameters) %in% random])
  )
  tape_random$simplify()
  tape_random$reorder()
  Hv = make_Hv(
    tape = tape_random,
    par = unlist(parameters[names(parameters) %in% random])
  )

  # Objective function
  nll = function(fixedvec){
    # Define fixed effects and assign to global environment
    env$par_vec[which(names(env$par_vec) %in% fixed)] = fixedvec
    tape_u$force.update()
    dfdu$force.update()

    if( is.null(env$inner_start) ){
      env$inner_start = env$par_vec[names(env$par_vec) %in% c(random,profile)]
      env$best = Inf
    }

    # Run inner and assign xhat to global environment
    out = optim(
      par = env$inner_start,
      fn = tape_u,
      gr = dfdu,
      method = "L-BFGS-B",
      control = list(trace=0, maxit = 1e3, factr = 1e-2)
    )

    # Define profiled effects, and assign to global environment
    which_profile = which( names(out$par) %in% profile )
    if( length(which_profile) > 0 ){
      env$par_vec[which(names(env$par_vec) %in% profile)] = out$par[which_profile]
      tape_random$force.update()
    }

    # Have to assign into env(Hv)$mle to evaluate at right point
    which_random = which( names(out$par) %in% random )
    if( length(which_random) > 0 ){
      attr(Hv,"env")$mle = out$par[which_random]
      env$logdet1_m = lanczos_logdet(
        Hv = Hv,
        k = k,
        m = m,
        n = length(which_random),
        seed = seed
      )
    }else{
      env$logdet1_m = rep(0, m)
    }

    # Get jnll
    jnll = tape_u( out$par )

    # Assemble laplace
    neglogmarglik = jnll - (0.5*length(out$par)*log(2*pi)) + 0.5*mean(env$logdet1_m)

    # Assign best and return
    if( neglogmarglik < env$best ){
      #env$inner_start = env$par_vec[names(env$par_vec) %in% c(random,profile)]
      env$inner_start = out$par
      env$best = neglogmarglik
    }
    return( neglogmarglik )
  }

  #
  out = list(
    par = unlist(parameters[fixed]),
    fn = nll,
    env = env,
    Hv = Hv
  )
  return(out)
}
