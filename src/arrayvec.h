#ifndef __ARRAYVEC_H__
#define __ARRAYVEC_H__

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

#include <strider.h>
using strider::make_strided;

template <size_t I>
using vec_type = array<double, I>;

template <size_t I>
using arrayvec = vector<vec_type<I>>;

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
  auto oi = back_inserter(*p);
  transform(begin(x), begin(x) + nr, oi,
            [&](const double& v)
            {
              vec_type<I> a;
              auto i = make_strided(&v, nr);
              copy(i, i + I, begin(a));
              return a;
            });
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
  transform(begin(*p), end(*p), begin(res), begin(res),
            [&](const array<double, I>& a, double& v)
            {
              auto i = make_strided(&v, p->size());
              copy(begin(a), end(a), i);
              return v;
            });
  return res;
}

template <size_t I>
NumericMatrix tuples_to_matrix_(List x, size_t a, size_t b)
{
  auto nr = b - a + 1;
  auto p = get_ptr<I>(x);
  if (b < a || p->size() < b + 1) stop("Invalid range");
  NumericMatrix res(nr, I);
  auto begin_ = begin(*p) + a,
       end_ = begin(*p) + b + 1;
  transform(begin_, end_, begin(res), begin(res),
            [&](const array<double, I>& u, double& v)
            {
              auto i = make_strided(&v, nr);
              copy(begin(u), end(u), i);
              return v;
            });
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

#endif // __ARRAYVEC_H__
