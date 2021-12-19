#include <Rcpp.h>
using Rcpp::as;
using Rcpp::stop;
using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::Function;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;

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

#ifndef NO_CXX17

#include "kdtools.h"
using namespace keittlab;
using kdtools::utils::median_part;
using kdtools::utils::iter_value_t;
using kdtools::utils::middle_of;
using kdtools::utils::n_best;

#include "strdist.h"
using namespace keittlab;

template<typename T>
int ncol(const T& x) {
  return x.length();
}

template<typename T>
int nrow(const T& x) {
  return Rf_xlength(SEXP(x[0]));
}

template<typename T, typename U>
bool not_in_range(const T& x, U n) {
  auto r = minmax_element(begin(x), end(x));
  return (*r.first < 1 || *r.second > n) ? true : false;
}

namespace {

std::string_view get_string(SEXP x, int i = 0) {
  return std::string_view(CHAR(STRING_ELT(x, i)));
}

Function Requal("=="), Rless("<"), Rdiff("-");

} // namespace

struct kd_less_df
{
  kd_less_df(const List& df, const IntegerVector& idx, int dim = 0, int count = 0)
    : m_df(df), m_idx(idx), m_dim(dim), m_ndim(m_idx.size()), m_count(count) {}

  kd_less_df next_dim(bool inc_count = false) const {
    return kd_less_df(m_df, m_idx,
                      (m_dim + 1) % m_ndim,
                      inc_count ? m_count + 1 : 0);
  }

  kd_less_df operator++() const { return next_dim(); }

  bool operator()(const int lhs, const int rhs) const {
    if (m_count == m_ndim) return false;
    auto col = SEXP(m_df[m_idx[m_dim] - 1]);
    switch(TYPEOF(col)) {
    case LGLSXP: {
      if (LOGICAL(col)[lhs] == LOGICAL(col)[rhs])
        return next_dim(true)(lhs, rhs);
      else
        return LOGICAL(col)[lhs] < LOGICAL(col)[rhs];
      break;
    }
    case REALSXP: {
      if (REAL(col)[lhs] == REAL(col)[rhs])
        return next_dim(true)(lhs, rhs);
      else
        return REAL(col)[lhs] < REAL(col)[rhs];
      break;
    }
    case INTSXP: {
      if (INTEGER(col)[lhs] == INTEGER(col)[rhs])
        return next_dim(true)(lhs, rhs);
      else
        return INTEGER(col)[lhs] < INTEGER(col)[rhs];
      break;
    }
    case STRSXP: {
      if (get_string(col, lhs) == get_string(col, rhs))
        return next_dim(true)(lhs, rhs);
      else
        return get_string(col, lhs) < get_string(col, rhs);
      break;
    }
    case VECSXP: {
      SEXP lhs_ = VECTOR_ELT(col, lhs),
        rhs_ = VECTOR_ELT(col, rhs);
      if (Requal(lhs_, rhs_))
        return next_dim(true)(lhs, rhs);
      else
        return Rless(lhs_, rhs_);
      break;
    }
    default: stop("Invalid column type");
    }
    return false;
  }
  const List& m_df;
  const IntegerVector& m_idx;
  int m_dim, m_ndim, m_count;
};

struct chck_nth_df
{
  chck_nth_df(const List& df, const IntegerVector& idx,
              const List& lower, const List& upper, int dim = 0)
    : m_df(df), m_lower(lower), m_upper(upper),
      m_idx(idx), m_dim(dim) {}

  chck_nth_df next_dim() const {
    return chck_nth_df(m_df, m_idx, m_lower, m_upper, (m_dim + 1) % m_idx.size());
  }

  chck_nth_df operator++() const { return next_dim(); }

  bool search_left(const int i) const {
    auto col = SEXP(m_df[m_idx[m_dim] - 1]),
      lwr = SEXP(m_lower[m_dim]);
    switch(TYPEOF(col)) {
    case LGLSXP: {
      if (!(LOGICAL(col)[i] < LOGICAL(lwr)[0])) return false;
      break;
    }
    case REALSXP: {
      if (!(REAL(col)[i] < REAL(lwr)[0])) return false;
      break;
    }
    case INTSXP: {
      if (!(INTEGER(col)[i] < INTEGER(lwr)[0])) return false;
      break;
    }
    case STRSXP: {
      if (!(get_string(col, i) < get_string(lwr))) return false;
      break;
    }
    case VECSXP: {
      SEXP lhs_ = VECTOR_ELT(col, i),
        rhs_ = VECTOR_ELT(lwr, 0);
      if (!Rless(lhs_, rhs_)) return false;
      break;
    }
    default: stop("Invalid column type");
    }
    return true;
  }

  bool search_right(const int i) const {
    auto col = SEXP(m_df[m_idx[m_dim] - 1]),
      upr = SEXP(m_upper[m_dim]);
    switch(TYPEOF(col)) {
    case LGLSXP: {
      if (LOGICAL(col)[i] < LOGICAL(upr)[0]) return true;
      break;
    }
    case REALSXP: {
      if (REAL(col)[i] < REAL(upr)[0]) return true;
      break;
    }
    case INTSXP: {
      if (INTEGER(col)[i] < INTEGER(upr)[0]) return true;
      break;
    }
    case STRSXP: {
      if (get_string(col, i) < get_string(upr)) return true;
      break;
    }
    case VECSXP: {
      SEXP lhs_ = VECTOR_ELT(col, i),
        rhs_ = VECTOR_ELT(upr, 0);
      if (Rless(lhs_, rhs_)) return true;
      break;
    }
    default: stop("Invalid column type");
    }
    return false;
  }

  const List& m_df, m_lower, m_upper;
  const IntegerVector& m_idx;
  int m_dim;
};

struct equal_nth_df
{
  equal_nth_df(const List& df, const IntegerVector& idx,
               const List& key, int dim = 0)
    : m_df(df), m_key(key), m_idx(idx), m_dim(dim) {}

  equal_nth_df next_dim() const {
    return equal_nth_df(m_df, m_idx, m_key, (m_dim + 1) % m_idx.size());
  }

  equal_nth_df operator++() const { return next_dim(); }

  bool operator()(const int i) const {
    auto col = SEXP(m_df[m_idx[m_dim] - 1]),
      key = SEXP(m_key[m_dim]);
    switch(TYPEOF(col)) {
    case LGLSXP: {
      if (LOGICAL(col)[i] == LOGICAL(key)[0]) return true;
      break;
    }
    case REALSXP: {
      if (REAL(col)[i] == REAL(key)[0]) return true;
      break;
    }
    case INTSXP: {
      if (INTEGER(col)[i] == INTEGER(key)[0]) return true;
      break;
    }
    case STRSXP: {
      if (get_string(col, i) == get_string(key)) return true;
      break;
    }
    case VECSXP: {
      SEXP lhs_ = VECTOR_ELT(col, i),
        rhs_ = VECTOR_ELT(key, 0);
      if (Requal(lhs_, rhs_)) return true;
      break;
    }
    default: stop("Invalid column type");
    }
    return false;
  }

  const List& m_df, m_key;
  const IntegerVector& m_idx;
  int m_dim;
};

struct dist_nth_df
{
  dist_nth_df(const List& df, const IntegerVector& idx,
              const NumericVector& w, const List& key, int dim = 0)
    : m_df(df), m_key(key), m_idx(idx), m_w(w), m_dim(dim) {}

  dist_nth_df next_dim() const {
    return dist_nth_df(m_df, m_idx, m_w, m_key, (m_dim + 1) % m_idx.size());
  }

  dist_nth_df operator++() const { return next_dim(); }

  double operator()(const int i) const {
    auto col = SEXP(m_df[m_idx[m_dim] - 1]),
      key = SEXP(m_key[m_dim]);
    switch(TYPEOF(col)) {
    case LGLSXP: {
      return LOGICAL(col)[i] == LOGICAL(key)[0] ? 0 : m_w[m_dim];
      break;
    }
    case REALSXP: {
      return m_w[m_dim] * std::abs(REAL(col)[i] - REAL(key)[0]);
      break;
    }
    case INTSXP: {
      if (Rf_inherits(col, "factor")) {
        return INTEGER(col)[i] == INTEGER(key)[0] ? 0 : m_w[m_dim];
      }
      return m_w[m_dim] * std::abs(INTEGER(col)[i] - INTEGER(key)[0]);
      break;
    }
    case STRSXP: {
      return m_w[m_dim] * strdist::levenshtein(get_string(col, i), get_string(key));
      break;
    }
    case VECSXP: {
      SEXP lhs_ = VECTOR_ELT(col, i), rhs_ = VECTOR_ELT(key, 0);
      return m_w[m_dim] * std::abs(as<double>(Rdiff(lhs_, rhs_)));
      break;
    }
    default: stop("Invalid column type");
    }
  }

  const List& m_df, m_key;
  const IntegerVector& m_idx;
  const NumericVector& m_w;
  int m_dim;
};

struct within_df {
  within_df(const List& df, const IntegerVector& idx,
            const List& lower, const List& upper)
    : m_df(df), m_lower(lower), m_upper(upper),
      m_idx(idx), m_ndim(m_idx.size()) {}

  bool operator()(const int i) const {
    for (int j = 0; j != m_ndim; ++j) {
      auto col = SEXP(m_df[m_idx[j] - 1]),
        l = SEXP(m_lower[j]), u = SEXP(m_upper[j]);
      switch(TYPEOF(col)) {
      case LGLSXP: {
        if (LOGICAL(col)[i] < LOGICAL(l)[0] ||
            LOGICAL(col)[i] >= LOGICAL(u)[0])
          return false;
        break;
      }
      case REALSXP: {
        if (REAL(col)[i] < REAL(l)[0] ||
            REAL(col)[i] >= REAL(u)[0])
          return false;
        break;
      }
      case INTSXP: {
        if (INTEGER(col)[i] < INTEGER(l)[0] ||
            INTEGER(col)[i] >= INTEGER(u)[0])
          return false;
        break;
      }
      case STRSXP: {
        if (get_string(col, i) < get_string(l) ||
            get_string(col, i) >= get_string(u))
          return false;
        break;
      }
      case VECSXP: {
        SEXP v = VECTOR_ELT(col, i);
        if (Rless(v, l) || (!Rless(u, v))) return false;
        break;
      }
      default: stop("Invalid column type");
      }
    }
    return true;
  }
  const List& m_df, m_lower, m_upper;
  const IntegerVector& m_idx;
  int m_ndim;
};

struct l2dist_df {
  l2dist_df(const List& df, const IntegerVector& idx, const NumericVector& w, const List& key)
    : m_df(df), m_key(key), m_idx(idx), m_w(w), m_ndim(m_idx.size()) {}

  double operator()(const int i) const {
    double ssq = 0;
    for (int j = 0; j != m_ndim; ++j) {
      auto col = SEXP(m_df[m_idx[j] - 1]), k = SEXP(m_key[j]);
      switch(TYPEOF(col)) {
      case LGLSXP: {
        ssq += LOGICAL(col)[i] == LOGICAL(k)[0] ? 0 : m_w[j];
        break;
      }
      case REALSXP: {
        ssq += m_w[j] * std::pow(REAL(col)[i] - REAL(k)[0], 2);
        break;
      }
      case INTSXP: {
        if (Rf_inherits(col, "factor")) {
          ssq += INTEGER(col)[i] == INTEGER(k)[0] ? 0 : m_w[j];
        } else {
          ssq += m_w[j] * std::pow(INTEGER(col)[i] - INTEGER(k)[0], 2);
        }
        break;
      }
      case STRSXP: {
        double d = strdist::levenshtein(get_string(col, i), get_string(k));
        ssq += m_w[j] * std::pow(d, 2);
        break;
      }
      case VECSXP: {
        SEXP lhs_ = VECTOR_ELT(col, i), rhs_ = VECTOR_ELT(k, 0);
        ssq += m_w[j] * std::pow(as<double>(Rdiff(lhs_, rhs_)), 2);
        break;
      }
      default: stop("Invalid column type");
      }
    }
    return std::sqrt(ssq);
  }
  const List& m_df, m_key;
  const IntegerVector& m_idx;
  const NumericVector& m_w;
  int m_ndim;
};

struct pdist_df {
  pdist_df(const List& df,
           const IntegerVector& idx,
           const NumericVector& w,
           const List& key, double p)
    : m_df(df), m_key(key), m_idx(idx), m_w(w), m_p(p), m_ndim(m_idx.size()) {}

  double operator()(const int i) const {
    double sum = 0;
    for (int j = 0; j != m_ndim; ++j) {
      auto col = SEXP(m_df[m_idx[j] - 1]), k = SEXP(m_key[j]);
      switch(TYPEOF(col)) {
      case LGLSXP: {
        sum += LOGICAL(col)[i] == LOGICAL(k)[0] ? 0 : m_w[j];
        break;
      }
      case REALSXP: {
        sum += m_w[j] * std::pow(std::abs(REAL(col)[i] - REAL(k)[0]), m_p);
        break;
      }
      case INTSXP: {
        if (Rf_inherits(col, "factor")) {
          sum += INTEGER(col)[i] == INTEGER(k)[0] ? 0 : m_w[j];
        } else {
          sum += m_w[j] * std::pow(std::abs(INTEGER(col)[i] - INTEGER(k)[0]), m_p);
        }
        break;
      }
      case STRSXP: {
        double d = strdist::levenshtein(get_string(col, i), get_string(k));
        sum += m_w[j] * std::pow(d, m_p);
        break;
      }
      case VECSXP: {
        SEXP lhs_ = VECTOR_ELT(col, i), rhs_ = VECTOR_ELT(k, 0);
        sum += m_w[j] * std::pow(std::abs(as<double>(Rdiff(lhs_, rhs_))), m_p);
        break;
      }
      default: stop("Invalid column type");
      }
    }
    return std::pow(sum, 1 / m_p);
  }
  const List& m_df, m_key;
  const IntegerVector& m_idx;
  const NumericVector& m_w;
  double m_p;
  int m_ndim;
};

bool type_mismatch(const List& df,
                   const IntegerVector& idx,
                   const List& lower,
                   const List& upper) {
  for (int i = 0; i != idx.size(); ++i) {
    int j = idx[i] - 1;
    auto c1 = SEXP(df[j]),
      c2 = SEXP(lower[i]),
      c3 = SEXP(upper[i]);
    if (TYPEOF(c1) != TYPEOF(c2) ||
        TYPEOF(c1) != TYPEOF(c3))
      return true;
  }
  return false;
}

bool type_mismatch(const List& df,
                   const IntegerVector& idx,
                   const List& value) {
  for (int i = 0; i != idx.size(); ++i) {
    int j = idx[i] - 1;
    auto c1 = SEXP(df[j]),
      c2 = SEXP(value[i]);
    if (TYPEOF(c1) != TYPEOF(c2))
      return true;
  }
  return false;
}

template <typename Iter, typename Pred>
void kd_order_df_(Iter first, Iter last, const Pred& pred)
{
  if (distance(first, last) > 1)
  {
    auto pivot = median_part(first, last, pred);
    kd_order_df_(next(pivot), last, ++pred);
    kd_order_df_(first, pivot, ++pred);
  }
}

template <typename Iter, typename Pred>
void kd_order_df_threaded(Iter first, Iter last, const Pred& pred,
                          int max_threads = std::thread::hardware_concurrency(),
                          int thread_depth = 1)
{
  if (distance(first, last) > 1)
  {
    auto pivot = median_part(first, last, pred);
    if ((1 << thread_depth) <= max_threads)
    {
      thread t(kd_order_df_threaded<Iter, Pred>,
               next(pivot), last, ++pred, max_threads, thread_depth + 1);
      kd_order_df_threaded<Iter, Pred>(first, pivot, ++pred, max_threads, thread_depth + 1);
      t.join();
    }
    else
    {
      kd_order_df_(next(pivot), last, ++pred);
      kd_order_df_(first, pivot, ++pred);
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
bool kd_is_sorted_df_(Iter first, Iter last, const Pred& pred)
{
  if (distance(first, last) < 2) return true;
  auto pivot = middle_of(first, last);
  return check_partition(first, pivot, last, pred) &&
    kd_is_sorted_df_(first, pivot, ++pred) &&
    kd_is_sorted_df_(next(pivot), last, ++pred);
}

template <typename Iter, typename Pred>
bool kd_is_sorted_df_threaded(Iter first, Iter last, const Pred& pred,
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
        res_left = kd_is_sorted_df_threaded(first, pivot, ++pred, max_threads, thread_depth + 1);
      });
      res_right = kd_is_sorted_df_threaded(next(pivot), last, ++pred, max_threads, thread_depth + 1);
      t.join();
      return res_left && res_right;
    }
    else
    {
      return kd_is_sorted_df_(first, pivot, ++pred) && kd_is_sorted_df_(next(pivot), last, ++pred);
    }
  } else {
    return false;
  }
}

template <typename Iter,
          typename OutIter,
          typename ChckNth,
          typename Within>
void kd_rq_df_(Iter first, Iter last, OutIter outp,
               const ChckNth& chck_nth,
               const Within& within)
{
  if (distance(first, last) > 32) {
    auto pivot = middle_of(first, last);
    if (within(*pivot)) *outp++ = *pivot;
    if (chck_nth.search_left(*pivot))
      kd_rq_df_(first, pivot, outp, ++chck_nth, within);
    if (chck_nth.search_right(*pivot))
      kd_rq_df_(next(pivot), last, outp, ++chck_nth, within);
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
          typename DistFun,
          typename QType>
void knn_(Iter first, Iter last,
          const EqualNth& equal_nth,
          const ChckNth& chck_nth,
          const DistNth& dist_nth,
          const DistFun& distf,
          QType& Q)
{
  switch(distance(first, last)) {
  case 1 : Q.add(distf(*first), first);
  case 0 : return; } // switch end
  auto pivot = middle_of(first, last);
  Q.add(distf(*pivot), pivot);
  if (equal_nth(*pivot)) {
    knn_(first, pivot, ++equal_nth, ++chck_nth, ++dist_nth, distf, Q);
    knn_(next(pivot), last, ++equal_nth, ++chck_nth, ++dist_nth, distf, Q);
  } else {
    auto search_left = !chck_nth.search_right(*pivot);
    if (search_left)
      knn_(first, pivot, ++equal_nth, ++chck_nth, ++dist_nth, distf, Q);
    else
      knn_(next(pivot), last, ++equal_nth, ++chck_nth, ++dist_nth, distf, Q);
    if (dist_nth(*pivot) <= Q.max_key())
    {
      if (search_left)
        knn_(next(pivot), last, ++equal_nth, ++chck_nth, ++dist_nth, distf, Q);
      else
        knn_(first, pivot, ++equal_nth, ++chck_nth, ++dist_nth, distf, Q);
    }
  }
}

template <typename Iter,
          typename EqualNth,
          typename ChckNth,
          typename DistNth,
          typename DistFun,
          typename QType>
void aknn_(Iter first, Iter last,
           const EqualNth& equal_nth,
           const ChckNth& chck_nth,
           const DistNth& dist_nth,
           const DistFun& distf,
           const double alpha,
           QType& Q)
{
  switch(distance(first, last)) {
  case 1 : Q.add(distf(*first), first);
  case 0 : return; } // switch end
  auto pivot = middle_of(first, last);
  Q.add(distf(*pivot), pivot);
  if (equal_nth(*pivot)) {
    aknn_(first, pivot, ++equal_nth, ++chck_nth, ++dist_nth, distf, alpha, Q);
    aknn_(next(pivot), last, ++equal_nth, ++chck_nth, ++dist_nth, distf, alpha, Q);
  } else {
    auto search_left = !chck_nth.search_right(*pivot);
    if (search_left)
      aknn_(first, pivot, ++equal_nth, ++chck_nth, ++dist_nth, distf, alpha, Q);
    else
      aknn_(next(pivot), last, ++equal_nth, ++chck_nth, ++dist_nth, distf, alpha, Q);
    if ((1 + alpha) * dist_nth(*pivot) <= Q.max_key())
    {
      if (search_left)
        aknn_(next(pivot), last, ++equal_nth, ++chck_nth, ++dist_nth, distf, alpha, Q);
      else
        aknn_(first, pivot, ++equal_nth, ++chck_nth, ++dist_nth, distf, alpha, Q);
    }
  }
}

#endif // NO_CXX17

// [[Rcpp::export]]
IntegerVector kd_order_df_no_validation(const List& df,
                                        const IntegerVector& idx,
                                        bool parallel = true) {
#ifdef NO_CXX17
  return R_NilValue;
#else
  IntegerVector x(nrow(df));
  iota(begin(x), end(x), 0);
  auto pred = kd_less_df(df, idx);
  if (parallel)
    kd_order_df_threaded(begin(x), end(x), pred);
  else
    kd_order_df_(begin(x), end(x), pred);
  return x + 1;
#endif
}

// [[Rcpp::export]]
IntegerVector kd_order_df(const List& df,
                          const IntegerVector& idx,
                          bool parallel = true) {
#ifdef NO_CXX17
  return R_NilValue;
#else
  if (ncol(df) < 1 || nrow(df) < 1)
    return IntegerVector();
  if (not_in_range(idx, ncol(df)))
    stop("Index out of range");
  return kd_order_df_no_validation(df, idx, parallel);
#endif
}

// [[Rcpp::export]]
bool kd_is_sorted_df_no_validation(const List& df,
                                   const IntegerVector& idx,
                                   bool parallel = true) {
#ifdef NO_CXX17
  return R_NilValue;
#else
  IntegerVector x(nrow(df));
  iota(begin(x), end(x), 0);
  auto pred = kd_less_df(df, idx);
  if (parallel)
    return kd_is_sorted_df_threaded(begin(x), end(x), pred);
  else
    return kd_is_sorted_df_(begin(x), end(x), pred);
#endif
}

// [[Rcpp::export]]
bool kd_is_sorted_df(const List& df,
                     const IntegerVector& idx,
                     bool parallel = true) {
#ifdef NO_CXX17
  return R_NilValue;
#else
  if (ncol(df) < 1 || nrow(df) < 1)
    stop("Invalid data frame");
  if (not_in_range(idx, ncol(df)))
    stop("Index out of range");
  return kd_is_sorted_df_no_validation(df, idx, parallel);
#endif
}

// [[Rcpp::export]]
std::vector<int> kd_rq_df_no_validation(const List& df,
                                        const IntegerVector& idx,
                                        const List& lower,
                                        const List& upper)
{
#ifdef NO_CXX17
  return std::vector<int>();
#else
  std::vector<int> x(nrow(df));
  iota(begin(x), end(x), 0);
  auto wi = within_df(df, idx, lower, upper);
  auto cn = chck_nth_df(df, idx, lower, upper);
  std::vector<int> res;
  auto oi = std::back_inserter(res);
  kd_rq_df_(begin(x), end(x), oi, cn, wi);
  for (auto& e : res) ++e;
  return res;
#endif
}

// [[Rcpp::export]]
std::vector<int> kd_rq_df(const List& df,
                          const IntegerVector& idx,
                          const List& lower,
                          const List& upper)
{
#ifdef NO_CXX17
  return std::vector<int>();
#else
  if (ncol(df) < 1 || nrow(df) < 1)
    stop("Empty data frame");
  if (not_in_range(idx, ncol(df)))
    stop("Index out of range");
  if (idx.size() != lower.size() ||
      idx.size() != upper.size())
    stop("Incorrect dimension of lower or upper bound");
  if (type_mismatch(df, idx, lower, upper))
    stop("Mismatched types in lower or upper bound");
  return kd_rq_df_no_validation(df, idx, lower, upper);
#endif
}

// [[Rcpp::export]]
List kd_nn_df_no_validation(const List& df,
                            const IntegerVector& idx,
                            const NumericVector& w,
                            const List& key,
                            const double a,
                            const double p,
                            const int n)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  auto m = std::min(n, nrow(df));
  std::vector<int> x(nrow(df));
  iota(begin(x), end(x), 0);
  auto equal_nth = equal_nth_df(df, idx, key);
  auto chck_nth = chck_nth_df(df, idx, key, key);
  auto distf = pdist_df(df, idx, w, key, p);
  auto dist_nth = dist_nth_df(df, idx, w, key);
  n_best<decltype(begin(x))> Q(m);
  if (a > 0)
    aknn_(begin(x), end(x), equal_nth, chck_nth, dist_nth, distf, a, Q);
  else
    knn_(begin(x), end(x), equal_nth, chck_nth, dist_nth, distf, Q);
  std::vector<std::pair<double, decltype(begin(x))>> out;
  auto oi = std::back_inserter(out);
  out.reserve(n);
  Q.copy_to(oi);
  IntegerVector loc(m);
  NumericVector dist(m);
  for (int i = 0; i != m; ++i) {
    loc[i] = distance(begin(x), out[i].second) + 1;
    dist[i] = out[i].first;
  }
  List res;
  res["index"] = loc;
  res["distance"] = dist;
  return res;
#endif
}

// [[Rcpp::export]]
List kd_nn_df(const List& df,
              const IntegerVector& idx,
              const NumericVector& w,
              const List& key,
              const double a,
              const double p,
              const int n)
{
#ifdef NO_CXX17
  return R_NilValue;
#else
  if (ncol(df) < 1 || nrow(df) < 1)
    stop("Empty data frame");
  if (not_in_range(idx, ncol(df)))
    stop("Index out of range");
  if (idx.size() != w.size())
    stop("Incorrect weights dimensions");
  if (idx.size() != key.size())
    stop("Incorrect dimension of key");
  if (type_mismatch(df, idx, key))
    stop("Mismatched types in key");
  return kd_nn_df_no_validation(df, idx, w, key, a, p, n);
#endif
}
