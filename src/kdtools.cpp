#include <Rcpp.h>
using Rcpp::NumericVector;
using Rcpp::XPtr;
using Rcpp::stop;
using Rcpp::List;
using Rcpp::as;

#include <iterator>
using std::begin;
using std::end;

#include <algorithm>
using std::copy;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <array>
using std::array;

#include "kdtools.h"

template <size_t I>
using arrayvec = vector<array<double, I>>;

template <typename T>
XPtr<T> make_xptr(T* x)
{
  return XPtr<T>(x);
}

template <size_t I, typename T>
XPtr<arrayvec<I>> get_ptr(const T& x)
{
  return as<XPtr<arrayvec<I>>>(x["xptr"]);
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
array<double, I> vec_to_array(const NumericVector& x)
{
  if (x.length() != I)
    stop("Invalid dimensions for value");
  array<double, I> y;
  copy(begin(x), end(x), begin(y));
  return y;
}

int arrayvec_dim(const List& x)
{
  if (!x.inherits("arrayvec"))
    stop("Expecting arrayvec object");
  return as<int>(x["ncol"]);
}

template <size_t I>
List kd_sort__(List x, bool inplace, bool parallel)
{
  auto p = get_ptr<I>(x);
  if (inplace)
  {
    if (parallel) kdtools::kd_sort_threaded(std::begin(*p), std::end(*p));
    else kdtools::kd_sort(std::begin(*p), std::end(*p));
    return x;
  }
  else
  {
    auto q = make_xptr(new arrayvec<I>(*p));
    if (parallel) kdtools::kd_sort_threaded(std::begin(*q), std::end(*q));
    else kdtools::kd_sort(std::begin(*q), std::end(*q));
    return wrap_ptr(q);
  }
}

// [[Rcpp::export]]
List kd_sort_(List x, bool inplace = false, bool parallel = false)
{
  switch(arrayvec_dim(x))
  {
  case 1: return kd_sort__<1>(x, inplace, parallel);
  case 2: return kd_sort__<2>(x, inplace, parallel);
  case 3: return kd_sort__<3>(x, inplace, parallel);
  case 4: return kd_sort__<4>(x, inplace, parallel);
  case 5: return kd_sort__<5>(x, inplace, parallel);
  case 6: return kd_sort__<6>(x, inplace, parallel);
  case 7: return kd_sort__<7>(x, inplace, parallel);
  case 8: return kd_sort__<8>(x, inplace, parallel);
  case 9: return kd_sort__<9>(x, inplace, parallel);
  default: stop("Invalid dimensions");
  }
}

template <size_t I>
List lex_sort__(List x, bool inplace)
{
  auto p = get_ptr<I>(x);
  if (inplace)
  {
    kdtools::lex_sort(std::begin(*p), std::end(*p));
    return x;
  }
  else
  {
    auto q = make_xptr(new arrayvec<I>(*p));
    kdtools::lex_sort(std::begin(*q), std::end(*q));
    return wrap_ptr(q);
  }
}

// [[Rcpp::export]]
List lex_sort_(List x, bool inplace = false)
{
  switch(arrayvec_dim(x))
  {
  case 1: return lex_sort__<1>(x, inplace);
  case 2: return lex_sort__<2>(x, inplace);
  case 3: return lex_sort__<3>(x, inplace);
  case 4: return lex_sort__<4>(x, inplace);
  case 5: return lex_sort__<5>(x, inplace);
  case 6: return lex_sort__<6>(x, inplace);
  case 7: return lex_sort__<7>(x, inplace);
  case 8: return lex_sort__<8>(x, inplace);
  case 9: return lex_sort__<9>(x, inplace);
  default: stop("Invalid dimensions");
  }
}

template <size_t I>
int kd_lower_bound__(List x, NumericVector v)
{
  auto p = get_ptr<I>(x);
  auto w = vec_to_array<I>(v);
  auto lv = kdtools::kd_lower_bound(begin(*p), end(*p), w);
  if (lv == end(*p)) return NA_INTEGER;
  return distance(begin(*p), lv) + 1;
}

// [[Rcpp::export]]
int kd_lower_bound_(List x, NumericVector value)
{
  switch(arrayvec_dim(x))
  {
  case 1: return kd_lower_bound__<1>(x, value);
  case 2: return kd_lower_bound__<2>(x, value);
  case 3: return kd_lower_bound__<3>(x, value);
  case 4: return kd_lower_bound__<4>(x, value);
  case 5: return kd_lower_bound__<5>(x, value);
  case 6: return kd_lower_bound__<6>(x, value);
  case 7: return kd_lower_bound__<7>(x, value);
  case 8: return kd_lower_bound__<8>(x, value);
  case 9: return kd_lower_bound__<9>(x, value);
  default: stop("Invalid dimensions");
  }
}

template <size_t I>
int kd_upper_bound__(List x, NumericVector v)
{
  auto p = get_ptr<I>(x);
  array<double, I> w;
  w = vec_to_array<I>(v);
  auto lv = kdtools::kd_upper_bound(begin(*p), end(*p), w);
  if (lv == end(*p)) return NA_INTEGER;
  return distance(begin(*p), lv) + 1;
}

// [[Rcpp::export]]
int kd_upper_bound_(List x, NumericVector value)
{
  switch(arrayvec_dim(x))
  {
  case 1: return kd_upper_bound__<1>(x, value);
  case 2: return kd_upper_bound__<2>(x, value);
  case 3: return kd_upper_bound__<3>(x, value);
  case 4: return kd_upper_bound__<4>(x, value);
  case 5: return kd_upper_bound__<5>(x, value);
  case 6: return kd_upper_bound__<6>(x, value);
  case 7: return kd_upper_bound__<7>(x, value);
  case 8: return kd_upper_bound__<8>(x, value);
  case 9: return kd_upper_bound__<9>(x, value);
  default: stop("Invalid dimensions");
  }
}

template <size_t I>
List kd_range_query__(List x, NumericVector lower, NumericVector upper)
{
  if (lower.length() != I || upper.length() != I)
    stop("Invalid dimensions for lower or upper");
  auto p = get_ptr<I>(x);
  auto q = make_xptr(new arrayvec<I>);
  auto oi = back_inserter(*q);
  auto l = vec_to_array<I>(lower),
    u = vec_to_array<I>(upper);
  kdtools::kd_range_query(begin(*p), end(*p), l, u, oi);
  return wrap_ptr(q);
}

// [[Rcpp::export]]
List kd_range_query_(List x, NumericVector lower, NumericVector upper)
{
  switch(arrayvec_dim(x))
  {
  case 1: return kd_range_query__<1>(x, lower, upper);
  case 2: return kd_range_query__<2>(x, lower, upper);
  case 3: return kd_range_query__<3>(x, lower, upper);
  case 4: return kd_range_query__<4>(x, lower, upper);
  case 5: return kd_range_query__<5>(x, lower, upper);
  case 6: return kd_range_query__<6>(x, lower, upper);
  case 7: return kd_range_query__<7>(x, lower, upper);
  case 8: return kd_range_query__<8>(x, lower, upper);
  case 9: return kd_range_query__<9>(x, lower, upper);
  default: stop("Invalid dimensions");
  }
}

template <size_t I>
int kd_nearest_neighbor__(List x, NumericVector v)
{
  auto p = get_ptr<I>(x);
  auto w = vec_to_array<I>(v);
  auto nn = kdtools::kd_nearest_neighbor(begin(*p), end(*p), w);
  if (nn >= end(*p)) stop("Search failed");
  return distance(begin(*p), nn) + 1;
}

// [[Rcpp::export]]
int kd_nearest_neighbor_(List x, NumericVector value)
{
  switch(arrayvec_dim(x))
  {
  case 1: return kd_nearest_neighbor__<1>(x, value);
  case 2: return kd_nearest_neighbor__<2>(x, value);
  case 3: return kd_nearest_neighbor__<3>(x, value);
  case 4: return kd_nearest_neighbor__<4>(x, value);
  case 5: return kd_nearest_neighbor__<5>(x, value);
  case 6: return kd_nearest_neighbor__<6>(x, value);
  case 7: return kd_nearest_neighbor__<7>(x, value);
  case 8: return kd_nearest_neighbor__<8>(x, value);
  case 9: return kd_nearest_neighbor__<9>(x, value);
  default: stop("Invalid dimensions");
  }
}

template <size_t I>
bool kd_binary_search__(List x, NumericVector v)
{
  auto p = get_ptr<I>(x);
  auto w = vec_to_array<I>(v);
  return kdtools::kd_binary_search(begin(*p), end(*p), w);
}

// [[Rcpp::export]]
bool kd_binary_search_(List x, NumericVector value)
{
  switch(arrayvec_dim(x))
  {
  case 1: return kd_binary_search__<1>(x, value);
  case 2: return kd_binary_search__<2>(x, value);
  case 3: return kd_binary_search__<3>(x, value);
  case 4: return kd_binary_search__<4>(x, value);
  case 5: return kd_binary_search__<5>(x, value);
  case 6: return kd_binary_search__<6>(x, value);
  case 7: return kd_binary_search__<7>(x, value);
  case 8: return kd_binary_search__<8>(x, value);
  case 9: return kd_binary_search__<9>(x, value);
  default: stop("Invalid dimensions");
  }
}

