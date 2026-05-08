
hutchinson_diag <-
function( Hq,
          m,
          n = length(attr(Hq,"env")$x0),
          seed = NULL ){

  if( !is.null(seed) ){
    set.seed(seed)
  }
  q_m = matrix( sample( c(-1,1), size = n*m, replace=TRUE), ncol = m )  # Rademacher vector
  diag_H = rowMeans(apply(q_m, MARGIN = 2, FUN = function(q) q * Hq(q) ))
  return(diag_H)
}

library(Matrix)
n = 1000
P = bandSparse( n = n, k = c(-1,1), diagonals = list(rep(0.5,n),rep(0.5,n)) )
Q = (Diagonal(n) - 0.9*t(P) ) %*% (Diagonal(n) - 0.9*P)
#L = array(rnorm(n^2), dim = c(n,n))
#Q = L %*% t(L)

Hq = function(x) (Q %*% x)[,1]

d1 = hutchinson_diag( Hq, m = 100, n = n )
d2 = diag(Q)
sqrt(mean((d1 - d2)^2)) / mean(abs(d2))
