#ifndef __INDEXED_H__
#define __INDEXED_H__

#include "arrayvec.h"

template <size_t I>
struct indexed
{
  static constexpr auto tuple_size = I;
  indexed(const array<double, I>& data, size_t index)
    : m_data(data), m_index(index) {}
  array<double, I> m_data;
  size_t m_index;
};

// WARNING Here be dragons!
namespace std {

template<size_t I, size_t N>
double& get(indexed<N> &x)
{
  return x.m_data[I];
}

template<size_t I, size_t N>
const double& get(const indexed<N> &x)
{
  return x.m_data[I];
}

template <size_t I>
class tuple_size<indexed<I>>
{
public:
  static constexpr auto value = indexed<I>::tuple_size;
};

}; // namespace std

template <size_t I>
using indexvec = vector<indexed<I>>;

template <size_t I, typename T>
XPtr<indexvec<I>> get_idx_ptr(const T& x)
{
  return as<XPtr<indexvec<I>>>(x["xptr"]);
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
  auto data_i = begin(*data_p);
  for (size_t i = 0; i != nr; ++i)
    indx_p->push_back(indexed<I>(*data_i++, i + 1));
  return wrap_idx_ptr(indx_p);
}

template <size_t I>
NumericMatrix indexed_to_matrix_(List x)
{
  auto p = get_idx_ptr<I>(x);
  NumericMatrix res(p->size(), I + 1);
  for (size_t i = 0; i != p->size(); ++i) {
    for (size_t j = 0; j != I; ++j)
      res(i, j) = (*p)[i].m_data[j];
    res(i, I) = (*p)[i].m_index;
  }
  return res;
}

template <size_t I>
NumericMatrix indexed_to_matrix_(List x, size_t a, size_t b)
{
  auto nr = b - a + 1;
  auto p = get_idx_ptr<I>(x);
  if (b < a || p->size() < b + 1) stop("Invalid range");
  NumericMatrix res(nr, I + 1);
  for (size_t i = a; i != b; ++i) {
    for (size_t j = 0; j != I; ++j)
      res(i, j) = (*p)[i].m_data[j];
    res(i, I) = (*p)[i].m_index;
  }
  return res;
}

inline
int indexvec_dim(const List& x)
{
  if (!x.inherits("indexvec"))
    stop("Expecting indexvec object");
  return as<int>(x["ncol"]);
}

#endif // __INDEXED_H__
