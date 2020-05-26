#include <Rcpp.h>
using Rcpp::as;
using Rcpp::stop;
using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::Function;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;

#include "kdtools.h"
using namespace keittlab;
using kdtools::utils::median_part;
using kdtools::utils::iter_value_t;
using kdtools::utils::middle_of;

#include <algorithm>
using std::end;
using std::next;
using std::swap;
using std::iota;
using std::begin;
using std::size_t;
using std::vector;
using std::thread;
using std::distance;
using std::minmax_element;
using std::partition_point;

// [[Rcpp::plugins(cpp17)]]

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

std::string_view get_string(SEXP x, int i) {
  return std::string_view(CHAR(STRING_ELT(x, i)));
}

Function Requal("=="), Rless("<");

struct kd_less_df
{
  kd_less_df(const List& df, const IntegerVector& idx, size_t dim = 0, size_t count = 0)
    : m_df(df), m_idx(idx), m_dim(dim), m_ndim(m_idx.size()), m_count(count) {}

  kd_less_df next_dim(bool inc_count = false) const {
    return kd_less_df(m_df, m_idx,
                      (m_dim + 1) % m_ndim,
                      inc_count ? m_count + 1 : 0);
  }

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
  size_t m_dim, m_ndim, m_count;
};

template <typename Iter, typename Pred>
void kd_order_df_(Iter first, Iter last, const Pred& pred)
{
  if (distance(first, last) > 1)
  {
    auto pivot = median_part(first, last, pred);
    kd_order_df_(next(pivot), last, pred.next_dim());
    kd_order_df_(first, pivot, pred.next_dim());
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
               next(pivot), last, pred.next_dim(), max_threads, thread_depth + 1);
      kd_order_df_threaded<Iter, Pred>(first, pivot, pred.next_dim(), max_threads, thread_depth + 1);
      t.join();
    }
    else
    {
      kd_order_df_(next(pivot), last, pred.next_dim());
      kd_order_df_(first, pivot, pred.next_dim());
    }
  }
}

// [[Rcpp::export]]
IntegerVector kd_order_df(const List& df,
                          const IntegerVector& idx,
                          bool parallel = true) {
  if (ncol(df) < 1 || nrow(df) < 1)
    return IntegerVector();
  if (not_in_range(idx, ncol(df)))
    stop("Index out of range");
  IntegerVector x(nrow(df));
  iota(begin(x), end(x), 0);
  auto pred = kd_less_df(df, idx);
  if (parallel)
    kd_order_df_threaded(begin(x), end(x), pred);
  else
    kd_order_df_(begin(x), end(x), pred);
  return x + 1;
}

struct less_nth_df
{
  less_nth_df(const List& df, const IntegerVector& idx,
              List lower, List upper, size_t dim = 0)
    : m_df(df), m_lower(lower), m_upper(upper),
      m_idx(idx), m_dim(dim) {}

  less_nth_df next_dim() const {
    return less_nth_df(m_df, m_idx, m_lower, m_upper, (m_dim + 1) % m_idx.size());
  }

  bool search_left(const int i) const {
    auto col = SEXP(m_df[m_idx[m_dim] - 1]),
      lwr = SEXP(m_lower[m_dim]);
    switch(TYPEOF(col)) {
    case LGLSXP: {
      if (LOGICAL(col)[i] < LOGICAL(lwr)[0]) return false;
      break;
    }
    case REALSXP: {
      if (REAL(col)[i] < REAL(lwr)[0]) return false;
      break;
    }
    case INTSXP: {
      if (INTEGER(col)[i] < INTEGER(lwr)[0]) return false;
      break;
    }
    case STRSXP: {
      if (get_string(col, i) < get_string(lwr, 0)) return false;
      break;
    }
    case VECSXP: {
      SEXP lhs_ = VECTOR_ELT(col, i),
        rhs_ = VECTOR_ELT(lwr, 0);
      if (Rless(lhs_, rhs_)) return false;
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
      if (get_string(col, i) < get_string(upr, 0)) return true;
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
  size_t m_dim;
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
        if (get_string(col, i) < get_string(l, 0) ||
            get_string(col, i) >= get_string(u, 0))
          return false;
        break;
      }
      default: stop("Invalid column type");
      }
    }
    return true;
  }
  const List& m_df, m_lower, m_upper;
  const IntegerVector& m_idx;
  size_t m_ndim;
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

template <typename Iter,
          typename OutIter,
          typename LessNth,
          typename Within,
          typename Pred>
void kd_rq_df_(Iter first, Iter last, OutIter outp,
               const Pred& pred,
               const LessNth& less_nth,
               const Within& within)
{
  if (distance(first, last) > 32) {
    auto pivot = middle_of(first, last);
    if (within(*pivot)) *outp++ = *pivot;
    if (less_nth.search_left(*pivot))
      kd_rq_df_(first, pivot, outp, pred.next_dim(), less_nth.next_dim(), within);
    if (less_nth.search_right(*pivot))
      kd_rq_df_(next(pivot), last, outp, pred.next_dim(), less_nth.next_dim(), within);
  } else {
    copy_if(first, last, outp, [&](const int x){
      return within(x);
    });
  }
  return;
}

// [[Rcpp::export]]
std::vector<int> kd_rq_df_no_validation(const List& df,
                                        const IntegerVector& idx,
                                        const List& lower,
                                        const List& upper)
{
  IntegerVector x(nrow(df));
  iota(begin(x), end(x), 0);
  auto pred = kd_less_df(df, idx);
  auto wi = within_df(df, idx, lower, upper);
  auto ln = less_nth_df(df, idx, lower, upper);
  std::vector<int> res;
  auto oi = std::back_inserter(res);
  kd_rq_df_(begin(x), end(x), oi, pred, ln, wi);
  std::transform(begin(res), end(res),
                 begin(res), [](int x){ return x + 1; });
  return res;
}

// [[Rcpp::export]]
std::vector<int> kd_rq_df(const List& df,
                          const IntegerVector& idx,
                          const List& lower,
                          const List& upper)
{
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
}

