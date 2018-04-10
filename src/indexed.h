#ifndef __INDEXED_H__
#define __INDEXED_H__

#include "kdtools.h"
using kdtools::key_value;

#include "arrayvec.h"

template <size_t I>
using indexvec = vector<key_value<array<double, I>, size_t>>;

template <size_t I, typename T>
XPtr<indexvec<I>> get_idx_ptr(const T& x)
{
  auto p = as<XPtr<indexvec<I>>>(x["xptr"]);
  if (!p) stop("Pointer is null");
  return p;
}

template <size_t I>
List wrap_idx_ptr(const XPtr<indexvec<I>>& q)
{
  List res;
  res["xptr"] = wrap(q);
  res["nrow"] = q->size();
  res["ncol"] = I + 1;
  res.attr("class") = "indexvec";
  return res;
}

template <size_t I>
List matrix_to_indexed_(const NumericMatrix& x)
{
  auto nr = x.nrow();
  auto data = matrix_to_tuples_<I>(x);
  auto data_p = get_av_ptr<I>(data);
  auto indx_p = make_xptr(new indexvec<I>);
  indx_p->reserve(nr);
  for (size_t i = 0; i != nr; ++i)
    indx_p->emplace_back(data_p->at(i), i + 1);
  return wrap_idx_ptr(indx_p);
}

template <size_t I>
NumericMatrix indexed_to_matrix_(List x)
{
  if (!x.inherits("indexvec"))
    stop("Expecting indexvec object");
  auto p = get_idx_ptr<I>(x);
  NumericMatrix res(p->size(), I + 1);
  for (size_t i = 0; i != p->size(); ++i) {
    for (size_t j = 0; j != I; ++j)
      res(i, j) = p->at(i).key().at(j);
    res(i, I) = p->at(i).value();
  }
  return res;
}

template <size_t I>
NumericMatrix indexed_to_matrix_(List x, size_t a, size_t b)
{
  if (!x.inherits("indexvec"))
    stop("Expecting indexvec object");
  auto nr = b - a + 1;
  auto p = get_idx_ptr<I>(x);
  if (b < a || p->size() < b + 1) stop("Invalid range");
  NumericMatrix res(nr, I + 1);
  for (size_t i = 0; i != nr; ++i) {
    for (size_t j = 0; j != I; ++j)
      res(i, j) = (p->at(i + a)).key().at(j);
    res(i, I) = (p->at(i + a)).value();
  }
  return res;
}

inline
int indexvec_dim(const List& x)
{
  return as<int>(x["ncol"]) - 1;
}

#endif // __INDEXED_H__
