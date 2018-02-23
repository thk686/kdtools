#include "arrayvec.h"
#include "kdtools.h"
using namespace kdtools;

template <size_t I>
List kd_sort__(List x, bool inplace, bool parallel)
{
  auto p = get_ptr<I>(x);
  if (inplace)
  {
    if (parallel) kd_sort_threaded(begin(*p), end(*p));
    else kd_sort(begin(*p), end(*p));
    return x;
  }
  else
  {
    auto q = make_xptr(new arrayvec<I>(*p));
    if (parallel) kd_sort_threaded(begin(*q), end(*q));
    else kd_sort(begin(*q), end(*q));
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
bool kd_is_sorted__(List x)
{
  auto p = get_ptr<I>(x);
  return kd_is_sorted(begin(*p), end(*p));
}

// [[Rcpp::export]]
bool kd_is_sorted_(List x)
{
  switch(arrayvec_dim(x))
  {
  case 1: return kd_is_sorted__<1>(x);
  case 2: return kd_is_sorted__<2>(x);
  case 3: return kd_is_sorted__<3>(x);
  case 4: return kd_is_sorted__<4>(x);
  case 5: return kd_is_sorted__<5>(x);
  case 6: return kd_is_sorted__<6>(x);
  case 7: return kd_is_sorted__<7>(x);
  case 8: return kd_is_sorted__<8>(x);
  case 9: return kd_is_sorted__<9>(x);
  default: stop("Invalid dimensions");
  }
}

template <size_t I>
List lex_sort__(List x, bool inplace)
{
  auto p = get_ptr<I>(x);
  if (inplace)
  {
    lex_sort(begin(*p), end(*p));
    return x;
  }
  else
  {
    auto q = make_xptr(new arrayvec<I>(*p));
    lex_sort(begin(*q), end(*q));
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
  auto lv = kd_lower_bound(begin(*p), end(*p), w);
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
  auto lv = kd_upper_bound(begin(*p), end(*p), w);
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
  auto p = get_ptr<I>(x);
  auto q = make_xptr(new arrayvec<I>);
  auto oi = back_inserter(*q);
  auto l = vec_to_array<I>(lower),
    u = vec_to_array<I>(upper);
  kd_range_query(begin(*p), end(*p), l, u, oi);
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
  auto nn = kd_nearest_neighbor(begin(*p), end(*p), w);
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
int kd_approx_nn__(List x, NumericVector v, double eps)
{
  auto p = get_ptr<I>(x);
  auto w = vec_to_array<I>(v);
  auto nn = kd_nearest_neighbor(begin(*p), end(*p), w, eps);
  if (nn >= end(*p)) stop("Search failed");
  return distance(begin(*p), nn) + 1;
}

// [[Rcpp::export]]
int kd_approx_nn_(List x, NumericVector value, double eps)
{
  switch(arrayvec_dim(x))
  {
  case 1: return kd_approx_nn__<1>(x, value, eps);
  case 2: return kd_approx_nn__<2>(x, value, eps);
  case 3: return kd_approx_nn__<3>(x, value, eps);
  case 4: return kd_approx_nn__<4>(x, value, eps);
  case 5: return kd_approx_nn__<5>(x, value, eps);
  case 6: return kd_approx_nn__<6>(x, value, eps);
  case 7: return kd_approx_nn__<7>(x, value, eps);
  case 8: return kd_approx_nn__<8>(x, value, eps);
  case 9: return kd_approx_nn__<9>(x, value, eps);
  default: stop("Invalid dimensions");
  }
}

template <size_t I>
bool kd_binary_search__(List x, NumericVector v)
{
  auto p = get_ptr<I>(x);
  auto w = vec_to_array<I>(v);
  return kd_binary_search(begin(*p), end(*p), w);
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

template <size_t I>
List kd_nearest_neighbors__(List x, NumericVector value, int n)
{
  auto p = get_ptr<I>(x);
  auto q = make_xptr(new arrayvec<I>);
  auto oi = back_inserter(*q);
  auto v = vec_to_array<I>(value);
  kd_nearest_neighbors(begin(*p), end(*p), v, n, oi);
  return wrap_ptr(q);
}

// [[Rcpp::export]]
List kd_nearest_neighbors_(List x, NumericVector value, int n)
{
  switch(arrayvec_dim(x))
  {
  case 1: return kd_nearest_neighbors__<1>(x, value, n);
  case 2: return kd_nearest_neighbors__<2>(x, value, n);
  case 3: return kd_nearest_neighbors__<3>(x, value, n);
  case 4: return kd_nearest_neighbors__<4>(x, value, n);
  case 5: return kd_nearest_neighbors__<5>(x, value, n);
  case 6: return kd_nearest_neighbors__<6>(x, value, n);
  case 7: return kd_nearest_neighbors__<7>(x, value, n);
  case 8: return kd_nearest_neighbors__<8>(x, value, n);
  case 9: return kd_nearest_neighbors__<9>(x, value, n);
  default: stop("Invalid dimensions");
  }
}
