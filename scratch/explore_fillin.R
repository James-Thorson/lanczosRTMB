

etree_from_chol <- function(L_factor) {
  n_col <- ncol(L_factor)
  parent <- integer(n_col)

  for (j in seq_len(n_col)) {
    col_start <- L_factor@p[j] + 1L
    col_end <- L_factor@p[j + 1L]
    if (col_start <= col_end) {
      rows <- L_factor@i[col_start:col_end] + 1L
      below_diag <- rows[rows > j]
      if (length(below_diag) > 0L) {
        parent[j] <- min(below_diag)
      }
    }
  }

  parent
}

summarize_cholesky_pattern <- function(Q) {
  chol_Q <- Cholesky(Q, perm = TRUE, LDL = FALSE)
  L <- as(chol_Q, "sparseMatrix")
  perm <- chol_Q@perm + 1L
  Q_perm <- as(Q[perm, perm], "CsparseMatrix")
  Q_pattern <- as(Q != 0, "CsparseMatrix")
  Q_perm_pattern <- as(Q_perm != 0, "CsparseMatrix")
  L_pattern <- as(L != 0, "CsparseMatrix")
  fill_pattern <- as((L_pattern + t(L_pattern)) != 0, "CsparseMatrix")
  fill_pattern <- as(fill_pattern & !Q_perm_pattern, "CsparseMatrix")

  list(
    dimension = dim(Q),
    nnz_Q = nnzero(Q_pattern),
    nnz_L = nnzero(L_pattern),
    nnz_fill = nnzero(fill_pattern),
    perm = perm,
    etree_parent = etree_from_chol(as(L, "dtCMatrix"))
  )
}

library(RTMB)
library(Matrix)
nx = 3
ndim = 3

Px = 0.9 * bandSparse( n = nx, k = c(1) )
Ix = Diagonal( n = nx )
Qx = (Ix - t(Px)) %*% (Ix - Px)

Q = Qx
for( i in seq_len(max(0,ndim-1)) ){
  Q = kronecker( Q, Qx )
}

summarize_cholesky_pattern(Q)

##########
# Incomplete
##########

X <- tril(Q)
## Store pattern of X
Q <- new("dgCMatrix", i=X@i, p=X@p, Dim=X@Dim, x=numeric(length(X@x)))
Q@x <- X@x
## Store forward pass Cholesky
L <- NULL
#logdet <- function(x) {
#  Q@x <- x
  L <<- .Call("tmb_ichol", t(Q), 1e-4, PACKAGE="TMB")
  L_pattern <- as(L != 0, "CsparseMatrix")
  Q_pattern <- as(Q != 0, "CsparseMatrix")
  fill_pattern <- as((L_pattern + t(L_pattern)) != 0, "CsparseMatrix")
  fill_pattern <- as(fill_pattern & !Q_pattern, "CsparseMatrix")
#  sum(log(diag(L)))
#}
