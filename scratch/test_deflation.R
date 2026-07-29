
# pak::pak( "James-Thorson/lanczosRTMB" )

library(Matrix)
library(lanczosRTMB)

nx = 100
q = sample( c(-1,1), size = nx, replace = TRUE )

rho = 0.99
P = bandSparse( n = nx, k = c(-1,1)) # , diagonals = list(rep(1,nx)) )
P = Diagonal( nx, 1/rowSums(P)) %*% P
I = Diagonal( nx )
H = (I - rho * t(P)) %*% (I - rho * P)

# True log-det
sum(log(eigen(H)$values))

###############
# SLQ on original HVP
###############

Hq = function(q,x) (H %*% q)[,1]

L = lanczos(
  Hq = Hq,
  q = q,
  k = 100,
  x = rep(0,nx),
  orthogonalize = TRUE # needed for nk = 100
)
Tri = tridiag(L$alpha, L$beta)

# Approximate log-det
sum(log(eigen(Tri)$values))

##############
# SLQ on deflated HVP
##############

V = L$Q %*% eigen(Tri)$vectors[,ncol(eigen(Tri)$vectors),drop=FALSE]
C = t(V) %*% H %*% V

#P = V %*% t(V)
# (I-P) %*% (I-P) = (I-P) by construction
# so compute qprime = (I-P) q
# w = H %*% vprime
# wprime = (I-P) w
# where wprime = (I-P) H (I-P)  %*% (I-P) q

L = lanczos(
  Hq = Hq,
  q = q,
  k = 100,
  V = V,
  x = rep(0,nx),
  orthogonalize = TRUE
)
Tri = tridiag(L$alpha, L$beta)

# Approximate log-det
sum(log(eigen(Tri)$values)) + log(C)
