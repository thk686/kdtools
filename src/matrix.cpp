#include <Rcpp.h>
using Rcpp::as;
using Rcpp::stop;
using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::Function;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::NumericMatrix;

#include <algorithm>
using std::end;
using std::next;
using std::swap;
using std::iota;
using std::begin;
using std::vector;
using std::distance;
using std::minmax_element;
using std::partition_point;

#include <thread>
using std::thread;

#include <cmath>

// #define USE_CIRCULAR_LEXICOGRAPHIC_COMPARE

#ifndef NO_CXX17

#include "kdtools.h"
using namespace keittlab;
using kdtools::utils::median_part;
using kdtools::utils::iter_value_t;
using kdtools::utils::middle_of;
using kdtools::utils::n_best;

#include "strdist.h"
using namespace keittlab;

template<typename T, typename U>
bool not_in_range(const T& x, U n) {
  auto r = minmax_element(begin(x), end(x));
  return (*r.first < 1 || *r.second > n) ? true : false;
}

namespace {

Function Requal("=="), Rless("<"), Rdiff("-");

};

#ifdef USE_CIRCULAR_LEXICOGRAPHIC_COMPARE

struct kd_less_mat
{
  kd_less_mat(const NumericMatrix& mat, const IntegerVector& idx, int dim = 0, int count = 0)
    : m_mat(mat), m_idx(idx), m_dim(dim), m_ndim(m_idx.size()), m_count(count) {}

  kd_less_mat next_dim(bool inc_count = false) const {
    return kd_less_mat(m_mat, m_idx,
                      (m_dim + 1) % m_ndim,
                      inc_count ? m_count + 1 : 0);
  }

  kd_less_mat operator++() const { return next_dim(); }

  bool operator()(const int lhs, const int rhs) const {
    if (m_count == m_ndim) return false;
    auto k = m_idx(m_dim) - 1;
    if (m_mat(lhs, k) == m_mat(rhs, k))
      return next_dim(true)(lhs, rhs);
    else
      return m_mat(lhs, k) < m_mat(rhs, k);
  }

  const NumericMatrix& m_mat;
  const IntegerVector& m_idx;
  int m_dim, m_ndim, m_count;
};

#else // (don't) USE_CIRCULAR_LEXICOGRAPHIC_COMPARE

struct kd_less_mat
{
  kd_less_mat(const NumericMatrix& mat, const IntegerVector& idx, int dim = 0)
    : m_mat(mat), m_idx(idx), m_dim(dim), m_ndim(m_idx.size()) {}

  kd_less_mat next_dim() const {
    return kd_less_mat(m_mat, m_idx, (m_dim + 1) % m_ndim);
  }

  kd_less_mat operator++() const {
    return next_dim();
  }

  bool operator()(const int lhs, const int rhs) const {
    auto k = m_idx(m_dim) - 1;
    return m_mat(lhs, k) < m_mat(rhs, k);
  }

  const NumericMatrix& m_mat;
  const IntegerVector& m_idx;
  int m_dim, m_ndim;
};

#endif // USE_CIRCULAR_LEXICOGRAPHIC_COMPARE

struct chck_nth_mat
{
  chck_nth_mat(const NumericMatrix& mat, const IntegerVector& idx,
              const NumericVector& lower, const NumericVector& upper, int dim = 0)
    : m_mat(mat), m_lower(lower), m_upper(upper),
      m_idx(idx), m_dim(dim) {}

  chck_nth_mat next_dim() const {
    return chck_nth_mat(m_mat, m_idx, m_lower, m_upper, (m_dim + 1) % m_idx.size());
  }

  chck_nth_mat operator++() const {
    return next_dim();
  }

  bool search_left(const int i) const {
    auto k = m_idx(m_dim) - 1;
    return(!(m_mat(i, k) < m_lower(k)));
  }

  bool search_right(const int i) const {
    auto k = m_idx(m_dim) - 1;
    return(m_mat(i, k) < m_upper(k));
  }

  const NumericMatrix& m_mat;
  const NumericVector& m_lower, m_upper;
  const IntegerVector& m_idx;
  int m_dim;
};

struct equal_nth_mat
{
  equal_nth_mat(const NumericMatrix& mat, const IntegerVector& idx,
                const NumericVector& key, int dim = 0)
    : m_mat(mat), m_key(key), m_idx(idx), m_dim(dim) {}

  equal_nth_mat next_dim() const {
    return equal_nth_mat(m_mat, m_idx, m_key, (m_dim + 1) % m_idx.size());
  }

  equal_nth_mat operator++() const { return next_dim(); }

  bool operator()(const int i) const {
    auto k = m_idx(m_dim) - 1;
    return(m_mat(i, k) == m_key(k));
  }

  const NumericMatrix& m_mat;
  const NumericVector& m_key;
  const IntegerVector& m_idx;
  int m_dim;
};

struct dist_nth_mat
{
  dist_nth_mat(const NumericMatrix& mat, const IntegerVector& idx,
               const NumericVector& key, int dim = 0)
    : m_mat(mat), m_key(key), m_idx(idx), m_dim(dim) {}

  dist_nth_mat next_dim() const {
    return dist_nth_mat(m_mat, m_idx, m_key, (m_dim + 1) % m_idx.size());
  }

  dist_nth_mat operator++() const { return next_dim(); }

  double operator()(const int i) const {
    auto k = m_idx(m_dim) - 1;
    return std::abs(m_mat(i, k) - m_key(k));
  }

  const NumericMatrix& m_mat;
  const NumericVector& m_key;
  const IntegerVector& m_idx;
  int m_dim;
};

struct within_mat {
  within_mat(const NumericMatrix& mat, const IntegerVector& idx,
            const NumericVector& lower, const NumericVector& upper)
    : m_mat(mat), m_lower(lower), m_upper(upper),
      m_idx(idx), m_ndim(m_idx.size()) {}

  bool operator()(const int i) const {
    for (int j = 0; j != m_ndim; ++j) {
      auto k = m_idx(j) - 1;
      if (m_mat(i, k) < m_lower(k) ||
          m_mat(i, k) >= m_upper(k))
        return false;
    }
    return true;
  }

  const NumericMatrix& m_mat;
  const NumericVector& m_lower, m_upper;
  const IntegerVector& m_idx;
  int m_ndim;
};

struct l2dist_mat {
  l2dist_mat(const NumericMatrix& mat, const IntegerVector& idx, const NumericVector& key)
    : m_mat(mat), m_key(key), m_idx(idx), m_ndim(m_idx.size()) {}

  double operator()(const int i) const {
    double ssq = 0;
    for (int j = 0; j != m_ndim; ++j) {
      auto k = m_idx(j) - 1;
      ssq += std::pow(m_mat(i, k) - m_key(k), 2);
    }
    return std::sqrt(ssq);
  }

  const NumericMatrix& m_mat;
  const NumericVector& m_key;
  const IntegerVector& m_idx;
  int m_ndim;
};

template <typename Iter, typename Pred>
void kd_order_mat_(Iter first, Iter last, const Pred& pred)
{
  if (distance(first, last) > 1)
  {
    auto pivot = median_part(first, last, pred);
    kd_order_mat_(next(pivot), last, ++pred);
    kd_order_mat_(first, pivot, ++pred);
  }
}

template <typename Iter, typename Pred>
void kd_order_mat_threaded(Iter first, Iter last, const Pred& pred,
                          int max_threads = std::thread::hardware_concurrency(),
                          int thread_depth = 1)
{
  if (distance(first, last) > 1)
  {
    auto pivot = median_part(first, last, pred);
    if ((1 << thread_depth) <= max_threads)
    {
      thread t(kd_order_mat_threaded<Iter, Pred>,
               next(pivot), last, ++pred, max_threads, thread_depth + 1);
      kd_order_mat_threaded<Iter, Pred>(first, pivot, ++pred, max_threads, thread_depth + 1);
      t.join();
    }
    else
    {
      kd_order_mat_(next(pivot), last, ++pred);
      kd_order_mat_(first, pivot, ++pred);
    }
  }
}

template <typename Iter, typename Pred>
bool check_partition(Iter first, Iter pivot, Iter last, Pred pred)
{
  while (first != pivot) if (pred(*pivot, *first++)) return false;
  while (first != last) if (pred(*first++, *pivot)) return false;
  return true;
}

template <typename Iter, typename Pred>
bool kd_is_sorted_mat_(Iter first, Iter last, const Pred& pred)
{
  if (distance(first, last) < 2) return true;
  auto pivot = middle_of(first, last);
  return check_partition(first, pivot, last, pred) &&
    kd_is_sorted_mat_(first, pivot, ++pred) &&
    kd_is_sorted_mat_(next(pivot), last, ++pred);
}

template <typename Iter, typename Pred>
bool kd_is_sorted_mat_threaded(Iter first, Iter last, const Pred& pred,
                               int max_threads = std::thread::hardware_concurrency(),
                               int thread_depth = 1)
{
  if (distance(first, last) < 2) return true;
  auto pivot = middle_of(first, last);
  if (check_partition(first, pivot, last, pred)) {
    if ((1 << thread_depth) <= max_threads)
    {
      bool res_left, res_right;
      thread t([=, &res_left](){
        res_left = kd_is_sorted_mat_threaded(first, pivot, ++pred, max_threads, thread_depth + 1);
      });
      res_right = kd_is_sorted_mat_threaded(next(pivot), last, ++pred, max_threads, thread_depth + 1);
      t.join();
      return res_left && res_right;
    }
    else
    {
      return kd_is_sorted_mat_(first, pivot, ++pred) && kd_is_sorted_mat_(next(pivot), last, ++pred);
    }
  } else {
    return false;
  }
}

template <typename Iter,
          typename OutIter,
          typename ChckNth,
          typename Within>
void kd_rq_mat_(Iter first, Iter last, OutIter outp,
               const ChckNth& chck_nth,
               const Within& within)
{
  if (distance(first, last) > 32) {
    auto pivot = middle_of(first, last);
    if (within(*pivot)) *outp++ = *pivot;
    if (chck_nth.search_left(*pivot))
      kd_rq_mat_(first, pivot, outp, ++chck_nth, within);
    if (chck_nth.search_right(*pivot))
      kd_rq_mat_(next(pivot), last, outp, ++chck_nth, within);
  } else {
    copy_if(first, last, outp, [&](const int x){
      return within(x);
    });
  }
  return;
}

template <typename Iter,
          typename EqualNth,
          typename ChckNth,
          typename DistNth,
          typename L2Dist,
          typename QType>
void knn_(Iter first, Iter last,
          const EqualNth& equal_nth,
          const ChckNth& chck_nth,
          const DistNth& dist_nth,
          const L2Dist& l2dist,
          QType& Q)
{
  switch(distance(first, last)) {
  case 1 : Q.add(l2dist(*first), first);
  case 0 : return; } // switch end
  auto pivot = middle_of(first, last);
  Q.add(l2dist(*pivot), pivot);
  if (equal_nth(*pivot)) {
    knn_(first, pivot, ++equal_nth, ++chck_nth, ++dist_nth, l2dist, Q);
    knn_(next(pivot), last, ++equal_nth, ++chck_nth, ++dist_nth, l2dist, Q);
  } else {
    auto search_left = !chck_nth.search_right(*pivot);
    if (search_left)
      knn_(first, pivot, ++equal_nth, ++chck_nth, ++dist_nth, l2dist, Q);
    else
      knn_(next(pivot), last, ++equal_nth, ++chck_nth, ++dist_nth, l2dist, Q);
    if (dist_nth(*pivot) <= Q.max_key())
    {
      if (search_left)
        knn_(next(pivot), last, ++equal_nth, ++chck_nth, ++dist_nth, l2dist, Q);
      else
        knn_(first, pivot, ++equal_nth, ++chck_nth, ++dist_nth, l2dist, Q);
    }
  }
}

#endif // NO_CXX17

// [[Rcpp::export]]
IntegerVector kd_order_mat_no_validation(const NumericMatrix& mat,
                                        const IntegerVector& idx,
                                        bool parallel = true) {
#ifndef NO_CXX17
  IntegerVector x(mat.nrow());
  iota(begin(x), end(x), 0);
  auto pred = kd_less_mat(mat, idx);
  if (parallel)
    kd_order_mat_threaded(begin(x), end(x), pred);
  else
    kd_order_mat_(begin(x), end(x), pred);
  return x + 1;
#else
  return IntegerVector();
#endif
}

// [[Rcpp::export]]
IntegerVector kd_order_mat(const NumericMatrix& mat,
                           const IntegerVector& idx,
                           bool parallel = true) {
#ifndef NO_CXX17
  if (mat.ncol() < 1 || mat.nrow() < 1)
    return IntegerVector();
  if (not_in_range(idx, mat.ncol()))
    stop("Index out of range");
  return kd_order_mat_no_validation(mat, idx, parallel);
#else
  return IntegerVector();
#endif
}

// [[Rcpp::export]]
bool kd_is_sorted_mat_no_validation(const NumericMatrix& mat,
                                    const IntegerVector& idx,
                                    bool parallel = true) {
#ifndef NO_CXX17
  IntegerVector x(mat.nrow());
  iota(begin(x), end(x), 0);
  auto pred = kd_less_mat(mat, idx);
  if (parallel)
    return kd_is_sorted_mat_threaded(begin(x), end(x), pred);
  else
    return kd_is_sorted_mat_(begin(x), end(x), pred);
#else
  return NA_LOGICAL;
#endif
}

// [[Rcpp::export]]
bool kd_is_sorted_mat(const NumericMatrix& mat,
                      const IntegerVector& idx,
                      bool parallel = true) {
#ifndef NO_CXX17
  if (mat.ncol() < 1 || mat.nrow() < 1)
    stop("Invalid input matrix");
  if (not_in_range(idx, mat.ncol()))
    stop("Index out of range");
  return kd_is_sorted_mat_no_validation(mat, idx, parallel);
#else
  return NA_LOGICAL;
#endif
}

// [[Rcpp::export]]
std::vector<int> kd_rq_mat_no_validation(const NumericMatrix& mat,
                                        const IntegerVector& idx,
                                        const NumericVector& lower,
                                        const NumericVector& upper)
{
#ifndef NO_CXX17
  std::vector<int> x(mat.nrow());
  iota(begin(x), end(x), 0);
  auto wi = within_mat(mat, idx, lower, upper);
  auto cn = chck_nth_mat(mat, idx, lower, upper);
  std::vector<int> res;
  auto oi = std::back_inserter(res);
  kd_rq_mat_(begin(x), end(x), oi, cn, wi);
  for (auto& e : res) ++e;
  return res;
#else
  return std::vector<int>();
#endif
}

// [[Rcpp::export]]
std::vector<int> kd_rq_mat(const NumericMatrix& mat,
                          const IntegerVector& idx,
                          const NumericVector& lower,
                          const NumericVector& upper)
{
#ifndef NO_CXX17
  if (mat.ncol() < 1 || mat.nrow() < 1)
    stop("Empty data frame");
  if (not_in_range(idx, mat.ncol()))
    stop("Index out of range");
  if (idx.size() != lower.size() ||
      idx.size() != upper.size())
    stop("Incorrect dimension of lower or upper bound");
  return kd_rq_mat_no_validation(mat, idx, lower, upper);
#else
  return std::vector<int>();
#endif
}

// [[Rcpp::export]]
std::vector<int> kd_nn_mat_no_validation(const NumericMatrix& mat,
                                        const IntegerVector& idx,
                                        const NumericVector& key,
                                        const int n)
{
#ifndef NO_CXX17
  std::vector<int> x(mat.nrow());
  iota(begin(x), end(x), 0);
  auto equal_nth = equal_nth_mat(mat, idx, key);
  auto chck_nth = chck_nth_mat(mat, idx, key, key);
  auto l2dist = l2dist_mat(mat, idx, key);
  auto dist_nth = dist_nth_mat(mat, idx, key);
  n_best<decltype(begin(x))> Q(n);
  knn_(begin(x), end(x), equal_nth, chck_nth, dist_nth, l2dist, Q);
  std::vector<int> res;
  auto oi = std::back_inserter(res);
  Q.copy_to(oi);
  for (auto& e : res) ++e;
  return res;
#else
  return std::vector<int>();
#endif
}

// [[Rcpp::export]]
std::vector<int> kd_nn_mat(const NumericMatrix& mat,
                           const IntegerVector& idx,
                           const NumericVector& key,
                           const int n)
{
#ifndef NO_CXX17
  if (mat.ncol() < 1 || mat.nrow() < 1)
    stop("Empty matrix");
  if (not_in_range(idx, mat.ncol()))
    stop("Index out of range");
  if (idx.size() != key.size())
    stop("Incorrect dimension of key");
  return kd_nn_mat_no_validation(mat, idx, key, n);
#else
  return std::vector<int>();
#endif
}

