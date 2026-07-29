
# pak::pak( "James-Thorson/lanczosRTMB" )

library(Matrix)
library(lanczosRTMB)

nx = 100
q = sample( c(-1,1), size = nx, replace = TRUE )
model = c( "SAR", "toy" )[2]

weight_logdet <- function(Tri, n_eff) {
  eT = eigen(Tri)
  w = eT$vectors[1,]^2
  n_eff * sum(w * log(eT$values))
}

if( model == "SAR" ){
  rho = 0.5
  P = bandSparse( n = nx, k = c(-1,1)) # , diagonals = list(rep(1,nx)) )
  P = Diagonal( nx, 1/rowSums(P)) %*% P
  I = Diagonal( nx )
  H = (I - rho * t(P)) %*% (I - rho * P)
}

if( model == "toy" ){
  # Random orthonormal basis (eigenvectors)
  Araw = matrix(rnorm(nx*nx), nx, nx)
  U = qr.Q(qr(Araw))   # nx x nx orthonormal

  # Spectrum: one tiny isolated eigenvalue, rest well-separated & fast-decaying
  lambda_bad = 1e-6                       # isolated small mode
  lambda_rest = 10^seq(2, 1, length.out = nx-1)   # geometric decay, e.g. 100 -> 1
  lambda = c(lambda_bad, lambda_rest)

  H = U %*% diag(lambda) %*% t(U)
  H = (H + t(H))/2   # symmetrize away roundoff
}

# True log-det
sum(log(eigen(H)$values))

###############
# SLQ on original HVP
###############
k = 10

Hq = function(q,x) (H %*% q)[,1]

L = lanczos(
  Hq = Hq,
  q = q,
  k = k,
  x = rep(0,nx),
  orthogonalize = TRUE # needed for nk = 100
)
Tri = tridiag(L$alpha, L$beta)

# Approximate log-det
weight_logdet(Tri, nx)

##############
# SLQ on deflated HVP
##############

extras = ncol(eigen(Tri)$vectors) - 0
V = L$Q %*% eigen(Tri)$vectors[,extras,drop=FALSE]
C = t(V) %*% H %*% V

#P = V %*% t(V)
# (I-P) %*% (I-P) = (I-P) by construction
# so compute qprime = (I-P) q
# w = H %*% vprime
# wprime = (I-P) w
# where wprime = (I-P) H (I-P)  %*% (I-P) q

L_deflated = lanczos(
  Hq = Hq,
  q = q,
  k = k,
  V = V,
  x = rep(0,nx),
  orthogonalize = TRUE
)
Tri_deflated = tridiag(L_deflated$alpha, L_deflated$beta)

# Approximate log-det
weight_logdet(Tri_deflated, nx - ncol(V)) + sum(log(eigen(C)$values))

##################
# Compare variance
##################

logdet_m = lanczos_logdet(
  Hq = Hq,
  k = k,
  m = 100,
  x = rep(0,nx),
  seed = 123,
  orthogonalize = TRUE,
  which_random = seq_len(nx)
)
var(logdet_m)

logdet_m = lanczos_logdet(
  Hq = Hq,
  k = k,
  m = 100,
  x = rep(0,nx),
  seed = 123,
  orthogonalize = TRUE,
  which_random = seq_len(nx),
  V = V
)
var(logdet_m)

