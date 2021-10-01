#include "arrayvec.h"

#ifndef NO_CXX17

#include "kdtools.h"
using namespace keittlab::kdtools;

template <size_t I>
List kd_sort__(List x, bool inplace, bool parallel)
{
  auto p = get_ptr<I>(x);
  if (inplace) {
    if (parallel) kd_sort_threaded(begin(*p), end(*p));
    else kd_sort(begin(*p), end(*p));
    return x;
  } else {
    auto q = make_xptr(new arrayvec<I>(*p));
    if (parallel) kd_sort_threaded(begin(*q), end(*q));
    else kd_sort(begin(*q), end(*q));
    return wrap_ptr(q);
  }
}

template <size_t I>
bool kd_is_sorted__(List x, bool parallel)
{
  auto p = get_ptr<I>(x);
  if (parallel)
    return kd_is_sorted_threaded(begin(*p), end(*p));
  else
    return kd_is_sorted(begin(*p), end(*p));
}

template <size_t I>
List lex_sort__(List x, bool inplace)
{
  auto p = get_ptr<I>(x);
  if (inplace) {
    lex_sort(begin(*p), end(*p));
    return x;
  } else {
    auto q = make_xptr(new arrayvec<I>(*p));
    lex_sort(begin(*q), end(*q));
    return wrap_ptr(q);
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

template <size_t I>
IntegerVector kd_rq_indices__(List x, NumericVector lower, NumericVector upper)
{
  auto p = get_ptr<I>(x);
  auto q = vector<av_iter<I>>();
  auto oi = back_inserter(q);
  auto l = vec_to_array<I>(lower),
    u = vec_to_array<I>(upper);
  kd_rq_iters(begin(*p), end(*p), l, u, oi);
  IntegerVector res(q.size());
  std::transform(begin(q), end(q), begin(res),
                 [&](av_iter<I> x){
                   return distance(p->begin(), x) + 1;
                 });
  return res;
}

template <size_t I>
List kd_rq_circular__(List x, NumericVector center, double radius)
{
  auto p = get_ptr<I>(x);
  auto q = make_xptr(new arrayvec<I>);
  auto oi = back_inserter(*q);
  auto cntr = vec_to_array<I>(center);
  kd_range_query(begin(*p), end(*p), cntr, radius, oi);
  return wrap_ptr(q);
}

template <size_t I>
IntegerVector kd_rqi_circular__(List x, NumericVector center, double radius)
{
  auto p = get_ptr<I>(x);
  auto q = vector<av_iter<I>>();
  auto oi = back_inserter(q);
  auto cntr = vec_to_array<I>(center);
  kd_rq_iters(begin(*p), end(*p), cntr, radius, oi);
  IntegerVector res(q.size());
  std::transform(begin(q), end(q), begin(res),
                 [&](av_iter<I> x){
                   return distance(p->begin(), x) + 1;
                 });
  return res;
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

template <size_t I>
bool kd_binary_search__(List x, NumericVector v)
{
  auto p = get_ptr<I>(x);
  auto w = vec_to_array<I>(v);
  return kd_binary_search(begin(*p), end(*p), w);
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

template <size_t I>
IntegerVector kd_nn_indices__(List x, NumericVector value, int n)
{
  auto p = get_ptr<I>(x);
  auto q = vector<av_iter<I>>();
  auto oi = back_inserter(q);
  auto v = vec_to_array<I>(value);
  kd_nn_iters(begin(*p), end(*p), v, n, oi);
  IntegerVector res(q.size());
  std::transform(begin(q), end(q), begin(res),
                 [&](av_iter<I> x){
                   return distance(p->begin(), x) + 1;
                 });
  return res;
}

template <size_t I>
List kd_nn_dist__(List x, NumericVector value, int n)
{
  auto p = get_ptr<I>(x);
  auto q = vector<std::pair<double, av_iter<I>>>();
  auto oi = back_inserter(q);
  auto v = vec_to_array<I>(value);
  q.reserve(n);
  kd_nn_dist(begin(*p), end(*p), v, n, oi);
  IntegerVector loc(n);
  NumericVector dist(n);
  for (int i = 0; i != n; ++i) {
    loc[i] = distance(p->begin(), q[i].second) + 1;
    dist[i] = q[i].first;
  }
  List res;
  res["index"] = loc;
  res["distance"] = dist;
  return res;
}

template <size_t I>
IntegerVector kd_order__(List x, bool inplace, bool parallel)
{
  auto p = get_ptr<I>(x);
  IntegerVector res(p->size());
  const vec_type<I>* p0 = &((*p)[0]);
  std::vector<vec_type<I>*> q(p->size());
  transform(begin(*p), end(*p), begin(q),
                 [](vec_type<I>& v){ return &v; });
  if (parallel) kd_sort_threaded(begin(q), end(q));
  else kd_sort(begin(q), end(q));
  transform(begin(q), end(q), begin(res),
            [p0](const vec_type<I>* v){
              return distance(p0, v) + 1;
            });
  if (inplace) {
    auto o = make_xptr(new arrayvec<I>);
    o->reserve(q.size());
    auto oi = back_inserter(*o);
    transform(begin(q), end(q), oi,
              [](const vec_type<I>* v){
                return *v;
              });
    x["xptr"] = wrap(o);
    p.release();
  }
  return res;
}

#endif // NO_CXX17

//' Check if package was compiled with circular comparisons
//'
//' @rdname utils
//' @export
// [[Rcpp::export]]
bool using_circular_lexicographical_compare() {
#ifdef USE_CIRCULAR_LEXICOGRAPHIC_COMPARE
  return true;
#else
  return false;
#endif
}

//' Check if C++ 17 was available when building package
//'
//' @rdname utils
//' @export
// [[Rcpp::export]]
bool has_cxx17() {
#ifdef NO_CXX17
  return false;
#else
  return true;
#endif
}

// [[Rcpp::export]]
IntegerVector kd_order_(List x, bool inplace, bool parallel)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  switch(arrayvec_dim(x)) {
  case 1: return kd_order__<1>(x, inplace, parallel);
  case 2: return kd_order__<2>(x, inplace, parallel);
  case 3: return kd_order__<3>(x, inplace, parallel);
  case 4: return kd_order__<4>(x, inplace, parallel);
  case 5: return kd_order__<5>(x, inplace, parallel);
  case 6: return kd_order__<6>(x, inplace, parallel);
  case 7: return kd_order__<7>(x, inplace, parallel);
  case 8: return kd_order__<8>(x, inplace, parallel);
  case 9: return kd_order__<9>(x, inplace, parallel);
  default: stop("Invalid dimensions");
  }
#endif
}

// [[Rcpp::export]]
bool kd_is_sorted_(List x, bool parallel)
{
#ifdef NO_CXX17
  return NA_LOGICAL;
#else
  switch(arrayvec_dim(x)) {
  case 1: return kd_is_sorted__<1>(x, parallel);
  case 2: return kd_is_sorted__<2>(x, parallel);
  case 3: return kd_is_sorted__<3>(x, parallel);
  case 4: return kd_is_sorted__<4>(x, parallel);
  case 5: return kd_is_sorted__<5>(x, parallel);
  case 6: return kd_is_sorted__<6>(x, parallel);
  case 7: return kd_is_sorted__<7>(x, parallel);
  case 8: return kd_is_sorted__<8>(x, parallel);
  case 9: return kd_is_sorted__<9>(x, parallel);
  default: stop("Invalid dimensions");
  }
#endif
}

// [[Rcpp::export]]
List kd_sort_(List x, bool inplace, bool parallel)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  switch(arrayvec_dim(x)) {
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
#endif
}

// [[Rcpp::export]]
List lex_sort_(List x, bool inplace)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  switch(arrayvec_dim(x)) {
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
#endif
}

// [[Rcpp::export]]
int kd_lower_bound_(List x, NumericVector value)
{
#ifdef NO_CXX17
  return NA_INTEGER;
#else
  switch(arrayvec_dim(x)) {
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
#endif
}

// [[Rcpp::export]]
int kd_upper_bound_(List x, NumericVector value)
{
#ifdef NO_CXX17
  return NA_INTEGER;
#else
  switch(arrayvec_dim(x)) {
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
#endif
}

// [[Rcpp::export]]
List kd_range_query_(List x, NumericVector lower, NumericVector upper)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  switch(arrayvec_dim(x)) {
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
#endif
}

// [[Rcpp::export]]
IntegerVector kd_rq_indices_(List x, NumericVector lower, NumericVector upper)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  switch(arrayvec_dim(x)) {
  case 1: return kd_rq_indices__<1>(x, lower, upper);
  case 2: return kd_rq_indices__<2>(x, lower, upper);
  case 3: return kd_rq_indices__<3>(x, lower, upper);
  case 4: return kd_rq_indices__<4>(x, lower, upper);
  case 5: return kd_rq_indices__<5>(x, lower, upper);
  case 6: return kd_rq_indices__<6>(x, lower, upper);
  case 7: return kd_rq_indices__<7>(x, lower, upper);
  case 8: return kd_rq_indices__<8>(x, lower, upper);
  case 9: return kd_rq_indices__<9>(x, lower, upper);
  default: stop("Invalid dimensions");
  }
#endif
}

// [[Rcpp::export]]
List kd_rq_circular_(List x, NumericVector center, double radius)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  switch(arrayvec_dim(x)) {
  case 1: return kd_rq_circular__<1>(x, center, radius);
  case 2: return kd_rq_circular__<2>(x, center, radius);
  case 3: return kd_rq_circular__<3>(x, center, radius);
  case 4: return kd_rq_circular__<4>(x, center, radius);
  case 5: return kd_rq_circular__<5>(x, center, radius);
  case 6: return kd_rq_circular__<6>(x, center, radius);
  case 7: return kd_rq_circular__<7>(x, center, radius);
  case 8: return kd_rq_circular__<8>(x, center, radius);
  case 9: return kd_rq_circular__<9>(x, center, radius);
  default: stop("Invalid dimensions");
  }
#endif
}

// [[Rcpp::export]]
IntegerVector kd_rqi_circular_(List x, NumericVector center, double radius)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  switch(arrayvec_dim(x)) {
  case 1: return kd_rqi_circular__<1>(x, center, radius);
  case 2: return kd_rqi_circular__<2>(x, center, radius);
  case 3: return kd_rqi_circular__<3>(x, center, radius);
  case 4: return kd_rqi_circular__<4>(x, center, radius);
  case 5: return kd_rqi_circular__<5>(x, center, radius);
  case 6: return kd_rqi_circular__<6>(x, center, radius);
  case 7: return kd_rqi_circular__<7>(x, center, radius);
  case 8: return kd_rqi_circular__<8>(x, center, radius);
  case 9: return kd_rqi_circular__<9>(x, center, radius);
  default: stop("Invalid dimensions");
  }
#endif
}

// [[Rcpp::export]]
int kd_nearest_neighbor_(List x, NumericVector value)
{
#ifdef NO_CXX17
  return NA_INTEGER;
#else
  switch(arrayvec_dim(x)) {
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
#endif
}

// [[Rcpp::export]]
bool kd_binary_search_(List x, NumericVector value)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  switch(arrayvec_dim(x)) {
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
#endif
}

// [[Rcpp::export]]
List kd_nearest_neighbors_(List x, NumericVector value, int n)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  switch(arrayvec_dim(x)) {
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
#endif
}

// [[Rcpp::export]]
IntegerVector kd_nn_indices_(List x, NumericVector value, int n)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  switch(arrayvec_dim(x)) {
  case 1: return kd_nn_indices__<1>(x, value, n);
  case 2: return kd_nn_indices__<2>(x, value, n);
  case 3: return kd_nn_indices__<3>(x, value, n);
  case 4: return kd_nn_indices__<4>(x, value, n);
  case 5: return kd_nn_indices__<5>(x, value, n);
  case 6: return kd_nn_indices__<6>(x, value, n);
  case 7: return kd_nn_indices__<7>(x, value, n);
  case 8: return kd_nn_indices__<8>(x, value, n);
  case 9: return kd_nn_indices__<9>(x, value, n);
  default: stop("Invalid dimensions");
  }
#endif
}

// [[Rcpp::export]]
List kd_nn_dist_(List x, NumericVector value, int n)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  switch(arrayvec_dim(x)) {
  case 1: return kd_nn_dist__<1>(x, value, n);
  case 2: return kd_nn_dist__<2>(x, value, n);
  case 3: return kd_nn_dist__<3>(x, value, n);
  case 4: return kd_nn_dist__<4>(x, value, n);
  case 5: return kd_nn_dist__<5>(x, value, n);
  case 6: return kd_nn_dist__<6>(x, value, n);
  case 7: return kd_nn_dist__<7>(x, value, n);
  case 8: return kd_nn_dist__<8>(x, value, n);
  case 9: return kd_nn_dist__<9>(x, value, n);
  default: stop("Invalid dimensions");
  }
#endif
}


