
namespace lanczosRTMB {

// Truncated conjugate gradient for solving y = A x for x given y and A;
// keeping z to easily add preconditioner later
template<class Type>
vector<Type> cg_solve(
      const vector<Type>& y,
      const Eigen::SparseMatrix<Type>& A,
      vector<Type> x,                           // initial guess
      int max_iter ){

  vector<Type> r = y - A * x;
  vector<Type> z = r;
  vector<Type> p = z;
  Type rz_old = r.dot(z);
  Type denom, rz_new, alpha, beta;
  for( int k = 0; k < max_iter; k++ ){
    vector<Type> Ap = A * p;
    alpha = rz_old / p.dot(Ap);
    x += alpha * p;
    r -= alpha * Ap;
    z = r;
    rz_new = r.dot(z);
    beta = rz_new / rz_old;
    p = z + beta * p;
    rz_old = rz_new;
  }
  return x;
}

}   // namespace lanczosRTMB

