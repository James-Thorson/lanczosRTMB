

#' @title
#' Approximate the log-density of a GMRF using Lanczos
#'
#' @description
#' Calculate log-density of a GMRF using Lanczos recursion
#'
#' @inheritParams lanczos
#' @inheritParams lanczos_logdet
#' @inheritParams RTMB::dgmrf
#'
#' @examples
#' library(RTMB)
#' library(Matrix)
#' library(lanczosRTMB)
#'
#' # Create a precision for a 1D AR1 process
#' nx = 1e5
#' P = 0.5 * bandSparse( n = nx, k = c(-1,1) )
#' Q = (Diagonal(nx) - 0.5*t(P) ) %*% (Diagonal(nx) - 0.5*P)
#' x = RTMB:::rgmrf0(n = 1, Q )[,1]
#'
#' # Exact log-density using RTMB
#' dgmrf(x, Q = Q, log = TRUE)
#'
#' # log-density using stochastic Lanczos quadrature
#' dgmrf_lanczos(x, Q = Q, k = 20, log = TRUE)
#'
#' # Log-density using incomplete Cholesky
#' dgmrf_icpl(x, Q = Q, log = TRUE)
#'
#' @export
dgmrf_lanczos <-
function( x,
          Q,
          k,
          mu = rep(0,length(x)),
          m = 3,
          log = FALSE,
          seed = NULL,
          orthogonalize = TRUE ){

  delta = x - mu
  Hq = function(x,...) (Q %*% x)[,1]
  env = list(
    which_random = seq_along(x),
    x = rep(0, length(x))
  )
  attr(Hq,"env") = env
  logdet_m = lanczos_logdet(
    Hq = Hq,
    k = k,
    m = m,
    orthogonalize = orthogonalize,
    seed = seed
  )

  nll = -0.5*sum(delta * (t(delta) %*% Q)) + 0.5*logdet_m - length(x)/2*log(2*pi)
  if(isTRUE(log)){return(nll)}else{return(exp(nll))}
}

#' @title
#' Approximate the log-determinant of a GMRF using incomplete Cholesky
#'
#' @export
logdet_icpl <- function(X) {
  ## Access only lower triangle
  X <- tril(X)
  ## Store pattern of X
  Q <- new("dgCMatrix", i=X@i, p=X@p, Dim=X@Dim, x=numeric(length(X@x)))
  ## Store forward pass Cholesky
  L <- NULL
  logdet <- function(x) {
    Q@x <- x
    L <<- .Call("tmb_ichol", t(Q), 1e-4, PACKAGE="TMB")
    sum(log(diag(L)))
  }
  dlogdet <- function(x, y, dy) {
    dL <- .Call("tmb_ldl_deriv", L, PACKAGE="TMB")
    i <- .Call("match_pattern", Q, dL, PACKAGE = "TMB")
    dL@x[i] * dy
  }
  RTMB::ADjoint(logdet, dlogdet)(X@x)
}

#' @title
#' Approximate the log-density of a GMRF using incomplete Cholesky
#'
#' @export
dgmrf_icpl = function( x, mu = 0, Q, k, v_m, log = FALSE ){
  logdetQ = logdet_icpl( X = Q )
  delta = x - mu
  logprior = -0.5*sum(delta * (t(delta) %*% Q)) + 0.5*logdetQ - length(x)/2*log(2*pi)
  if(isTRUE(log)){return(logprior)}else{return(exp(logprior))}
}




