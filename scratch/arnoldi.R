arnoldi <-
function( Hq,
          q,
          k,
          x = attr(Hq,"env")$x0,
          orthogonalize = TRUE,
          tol = 1e-12 ) {

  n = length(q)

  # Arnoldi basis
  Q = matrix(0, n, k + 1)

  # Upper Hessenberg matrix
  H = matrix(0, k + 1, k)

  # Normalize initial vector
  q = q / sqrt(sum(q*q))
  Q[,1] = q

  m = k   # actual number of Arnoldi steps performed

  for (j in 1:k) {

    # Apply operator
    w = Hq(Q[,j], x)

    # Modified Gram-Schmidt
    for (i in 1:j) {
      H[i,j] = sum(Q[,i] * w)
      w = w - H[i,j] * Q[,i]
    }

    if( isTRUE(orthogonalize) ){
      # Optional second pass for stability
      for (i in 1:j) {
        h_corr = sum(Q[,i] * w)
        H[i,j] = H[i,j] + h_corr
        w = w - h_corr * Q[,i]
      }
    }

    # Subdiagonal Hessenberg entry
    H[j+1,j] = sqrt(sum(w*w))

    # Check convergence
    if(H[j+1,j] < tol) {
      m = j
      break
    }

    # Normalize next basis vector
    Q[,j+1] = w / H[j+1,j]
  }

  list(
    Q = Q[, seq_len(m + 1), drop = FALSE],
    H = H[seq_len(m + 1), seq_len(m), drop = FALSE],
    m = m
  )
}
