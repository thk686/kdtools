#include "arrayvec.h"

//' Convert a matrix to a vector of arrays
//'
//' @param x object to be converted
//'
//' @details The algorithms in kdtools can accept either matrices or an
//' \link{arrayvec} object. When a matrix is passed, it is converted to
//' an arrayvec object internally and the results are converted back to
//' a matrix. For optimal performance, pre-convert matrices.
//'
//' @examples
//' if (has_cxx17()) {
//' x = matrix(1:10, 5)
//' y = matrix_to_tuples(x)
//' str(x)
//' str(y)
//' y[1:2, ]
//' }
//'
//' @rdname convert
//' @export
// [[Rcpp::export]]
List matrix_to_tuples(const NumericMatrix& x)
{
#ifdef HAS_CXX17
  switch(x.ncol())
  {
  case 1: return matrix_to_tuples_<1>(x);
  case 2: return matrix_to_tuples_<2>(x);
  case 3: return matrix_to_tuples_<3>(x);
  case 4: return matrix_to_tuples_<4>(x);
  case 5: return matrix_to_tuples_<5>(x);
  case 6: return matrix_to_tuples_<6>(x);
  case 7: return matrix_to_tuples_<7>(x);
  case 8: return matrix_to_tuples_<8>(x);
  case 9: return matrix_to_tuples_<9>(x);
  default: stop("Invalid dimensions");
  }
#else
  return List();
#endif
}

//' @rdname convert
//' @export
// [[Rcpp::export]]
NumericMatrix tuples_to_matrix(List x)
{
#ifdef HAS_CXX17
  if (!x.inherits("arrayvec"))
    stop("Expecting arrayvec object");
  switch(arrayvec_dim(x))
  {
  case 1: return tuples_to_matrix_<1>(x);
  case 2: return tuples_to_matrix_<2>(x);
  case 3: return tuples_to_matrix_<3>(x);
  case 4: return tuples_to_matrix_<4>(x);
  case 5: return tuples_to_matrix_<5>(x);
  case 6: return tuples_to_matrix_<6>(x);
  case 7: return tuples_to_matrix_<7>(x);
  case 8: return tuples_to_matrix_<8>(x);
  case 9: return tuples_to_matrix_<9>(x);
  default: stop("Invalid dimensions");
  }
#else
  return NumericMatrix();
#endif
}

// [[Rcpp::export]]
NumericMatrix tuples_to_matrix_rows(List x, int a, int b)
{
#ifdef HAS_CXX17
  if (!x.inherits("arrayvec"))
    stop("Expecting arrayvec object");
  switch(arrayvec_dim(x))
  {
  case 1: return tuples_to_matrix_<1>(x, a, b);
  case 2: return tuples_to_matrix_<2>(x, a, b);
  case 3: return tuples_to_matrix_<3>(x, a, b);
  case 4: return tuples_to_matrix_<4>(x, a, b);
  case 5: return tuples_to_matrix_<5>(x, a, b);
  case 6: return tuples_to_matrix_<6>(x, a, b);
  case 7: return tuples_to_matrix_<7>(x, a, b);
  case 8: return tuples_to_matrix_<8>(x, a, b);
  case 9: return tuples_to_matrix_<9>(x, a, b);
  default: stop("Invalid dimensions");
  }
#else
  return NumericMatrix();
#endif
}

