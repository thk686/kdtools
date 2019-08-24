#include <Rcpp.h>
using Rcpp::as;
using Rcpp::stop;
using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::Function;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;

#include "kdtools.h"
using kdtools::utils::median_part;

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

std::string_view get_string(SEXP x, int i) {
  return std::string_view(CHAR(STRING_ELT(x, i)));
}

Function Requal("=="), Rless("<");

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

template <typename Iter, typename Pred>
void kd_order_df_threaded(Iter first, Iter last, Pred pred,
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
IntegerVector kd_order_df(List df, IntegerVector idx, bool parallel = true)
{
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

/*
template <typename Iter,  typename Pred>
Iter find_pivot(Iter first, Iter last, Pred pred)
{
  using T = iter_value_t<Iter>;
  auto pivot = middle_of(first, last);
  return partition_point(first, pivot, [&](const T& x){
    return pred(x, *pivot);
  });
}

template <typename Iter,
          typename TupleType,
          typename OutIter,
          typename Pred>
void kd_range_query(Iter first, Iter last,
                    const TupleType& lower,
                    const TupleType& upper,
                    OutIter outp, Pred pred)
{
  if (distance(first, last) > 32) {
    auto pivot = find_pivot(first, last, pred);
    if (within(*pivot, lower, upper)) *outp++ = *pivot;
    if (!pred(*pivot, lower)) // search left
      kd_range_query(first, pivot, lower, upper, outp, pred);
    if (pred(*pivot, upper)) // search right
      kd_range_query(next(pivot), last, lower, upper, outp, pred);
  } else {
    copy_if(first, last, outp, [&](const TupleType& x){
      return within(x, lower, upper);
    });
  }
  return;
}
*/


