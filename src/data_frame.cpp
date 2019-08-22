#include <Rcpp.h>
using Rcpp::as;
using Rcpp::stop;
using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;

#include "kdtools.h"
using kdtools::detail::median_part;

#include <algorithm>
using std::end;
using std::next;
using std::swap;
using std::iota;
using std::begin;
using std::size_t;
using std::vector;
using std::distance;
using std::minmax_element;

template<typename T>
int ncol(const T& x) {
  return x.length();
}

template<typename T>
int nrow(T& x) {
  return Rf_xlength(x[0].get());
}

template<typename T, typename U>
bool not_in_range(const T& x, U n) {
  auto r = minmax_element(begin(x), end(x));
  return (*r.first < 1 || *r.second > n) ? true : false;
}

struct kd_less_df
{
  kd_less_df(List& df, IntegerVector& idx, size_t dim = 0, size_t count = 0)
    : m_df(df), m_idx(idx), m_dim(dim), m_ndim(m_idx.size()), m_count(count) {}
  kd_less_df next_dim(bool inc_count = false) {
    return kd_less_df(m_df, m_idx,
                      (m_dim + 1) % m_ndim,
                      inc_count ? m_count + 1 : 0);
  }
  bool operator()(const int lhs, const int rhs)
  {
    if (m_count == m_ndim) return false;
    auto col = m_df[m_idx[m_dim] - 1].get();
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
    default: stop("Invalid column type");
    }
    return false;
  }
  List& m_df;
  IntegerVector& m_idx;
  size_t m_dim, m_ndim, m_count;
};

template <typename Iter, typename Pred>
void kd_order_df_(Iter first, Iter last, Pred pred)
{
  if (distance(first, last) > 1)
  {
    auto pivot = median_part(first, last, pred);
    kd_order_df_(next(pivot), last, pred.next_dim());
    kd_order_df_(first, pivot, pred.next_dim());
  }
}

//' @export
// [[Rcpp::export]]
IntegerVector kd_order_df(List df, IntegerVector idx)
{
  if (ncol(df) < 1 || nrow(df) < 1)
    return IntegerVector();
  if (not_in_range(idx, ncol(df)))
    stop("Index out of range");
  IntegerVector x(nrow(df));
  iota(begin(x), end(x), 0);
  auto pred = kd_less_df(df, idx);
  kd_order_df_(begin(x), end(x), pred);
  return x + 1;
}
