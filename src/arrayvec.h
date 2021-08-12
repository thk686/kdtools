#ifndef ARRAYVEC_H
#define ARRAYVEC_H

#include <Rcpp.h>
using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::NumericMatrix;
using Rcpp::stop;
using Rcpp::XPtr;
using Rcpp::List;
using Rcpp::wrap;
using Rcpp::as;

#include <algorithm>
using std::transform;
using std::copy;

#include <array>
using std::array;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <iterator>
using std::back_inserter;
using std::distance;
using std::begin;
using std::end;

#ifndef NO_CXX17

template <size_t I>
using vec_type = array<double, I>;

template <size_t I>
using arrayvec = vector<vec_type<I>>;

template <size_t I>
using av_iter = typename arrayvec<I>::iterator;

template <typename T>
XPtr<T> make_xptr(T* x)
{
  return XPtr<T>(x);
}

template <size_t I>
List wrap_ptr(const XPtr<arrayvec<I>>& q)
{
  List res;
  res["xptr"] = wrap(q);
  res["nrow"] = q->size();
  res["ncol"] = I;
  res.attr("class") = "arrayvec";
  return res;
}

template <size_t I>
List matrix_to_tuples_(const NumericMatrix& x)
{
  auto nr = x.nrow();
  auto p = make_xptr(new arrayvec<I>);
  p->reserve(nr);
  for (auto i = 0; i != nr; ++i) {
    vec_type<I> a;
    for (auto j = 0; j != I; ++j) a[j] = x(i, j);
    p->push_back(a);
  }
  return wrap_ptr(p);
}

template <size_t I, typename T>
XPtr<arrayvec<I>> get_ptr(const T& x)
{
  return as<XPtr<arrayvec<I>>>(x["xptr"]);
}

template <size_t I>
NumericMatrix tuples_to_matrix_(List x)
{
  auto p = get_ptr<I>(x);
  NumericMatrix res(p->size(), I);
  for (auto i = 0; i != res.nrow(); ++i)
    for (auto j = 0; j != I; ++j)
      res(i, j) = (*p)[i][j];
  return res;
}

template <size_t I>
NumericMatrix tuples_to_matrix_(List x, size_t a, size_t b)
{
  auto nr = b - a + 1;
  auto p = get_ptr<I>(x);
  if (a < 1 || b < a || p->size() < b) stop("Invalid range");
  NumericMatrix res(nr, I);
  for (auto i = a; i != b + 1; ++i)
    for (auto j = 0; j != I; ++j)
      res(i - 1, j - 1) = (*p)[i - 1][j - 1];
  return res;
}

template <size_t I>
vec_type<I> vec_to_array(const NumericVector& x)
{
  if (x.length() != I)
    stop("Invalid dimensions for value");
  vec_type<I> y;
  copy(begin(x), end(x), begin(y));
  return y;
}

inline
int arrayvec_dim(const List& x)
{
  if (!x.inherits("arrayvec"))
    stop("Expecting arrayvec object");
  return as<int>(x["ncol"]);
}

#endif // NO_CXX17

#endif // ARRAYVEC_H
