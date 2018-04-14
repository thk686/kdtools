#include "indexed.h"

//' Convert a matrix to a vector of arrays
//'
//' @examples
//' x = matrix(1:10, 5)
//' y = matrix_to_indexed(x)
//' str(x)
//' str(y)
//' y[1:2, ]
//'
//' @rdname convert
//' @export
// [[Rcpp::export]]
List matrix_to_indexed(const NumericMatrix& x)
{
  switch(x.ncol())
  {
  case 1: return matrix_to_indexed_<1>(x);
  case 2: return matrix_to_indexed_<2>(x);
  case 3: return matrix_to_indexed_<3>(x);
  case 4: return matrix_to_indexed_<4>(x);
  case 5: return matrix_to_indexed_<5>(x);
  case 6: return matrix_to_indexed_<6>(x);
  case 7: return matrix_to_indexed_<7>(x);
  case 8: return matrix_to_indexed_<8>(x);
  case 9: return matrix_to_indexed_<9>(x);
  default: stop("Invalid dimensions");
  }
}

//' @rdname convert
//' @export
// [[Rcpp::export]]
NumericMatrix indexed_to_matrix(List x)
{
  if (!x.inherits("indexvec"))
    stop("Expecting indexvec object");
  switch(indexvec_dim(x))
  {
  case 1: return indexed_to_matrix_<1>(x);
  case 2: return indexed_to_matrix_<2>(x);
  case 3: return indexed_to_matrix_<3>(x);
  case 4: return indexed_to_matrix_<4>(x);
  case 5: return indexed_to_matrix_<5>(x);
  case 6: return indexed_to_matrix_<6>(x);
  case 7: return indexed_to_matrix_<7>(x);
  case 8: return indexed_to_matrix_<8>(x);
  case 9: return indexed_to_matrix_<9>(x);
  default: stop("Invalid dimensions");
  }
}

// [[Rcpp::export]]
NumericMatrix indexed_to_matrix_rows(List x, int a, int b)
{
  if (!x.inherits("indexvec"))
    stop("Expecting indexvec object");
  switch(indexvec_dim(x))
  {
  case 1: return indexed_to_matrix_<1>(x, a, b);
  case 2: return indexed_to_matrix_<2>(x, a, b);
  case 3: return indexed_to_matrix_<3>(x, a, b);
  case 4: return indexed_to_matrix_<4>(x, a, b);
  case 5: return indexed_to_matrix_<5>(x, a, b);
  case 6: return indexed_to_matrix_<6>(x, a, b);
  case 7: return indexed_to_matrix_<7>(x, a, b);
  case 8: return indexed_to_matrix_<8>(x, a, b);
  case 9: return indexed_to_matrix_<9>(x, a, b);
  default: stop("Invalid dimensions");
  }
}
