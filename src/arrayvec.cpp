#include <Rcpp.h>
using Rcpp::NumericMatrix;
using Rcpp::stop;
using Rcpp::XPtr;
using Rcpp::List;
using Rcpp::wrap;
using Rcpp::as;

#include <algorithm>
using std::transform;
using std::copy;

#include <array>
using std::array;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <iterator>
using std::back_inserter;
using std::begin;
using std::end;

#include <strider.h>
using strider::make_strided;

template <size_t I>
using arrayvec = vector<array<double, I>>;

template <typename T>
XPtr<T> make_xptr(T* x)
{
  return XPtr<T>(x);
}

template <size_t I>
List matrix_to_tuples_(const NumericMatrix& x)
{
  auto nr = x.nrow();
  auto p = make_xptr(new arrayvec<I>);
  p->reserve(nr);
  auto oi = back_inserter(*p);
  transform(begin(x), begin(x) + nr, oi,
            [&](const double& v)
            {
              array<double, I> a;
              auto i = make_strided(&v, nr);
              copy(i, i + I, begin(a));
              return a;
            });
  List res;
  res["xptr"] = wrap(p);
  res["nrow"] = wrap(nr);
  res["ncol"] = wrap(I);
  res.attr("class") = "arrayvec";
  return res;
}

template <size_t I, typename T>
XPtr<arrayvec<I>> get_ptr(const T& x)
{
  return as<XPtr<arrayvec<I>>>(x["xptr"]);
}

template <size_t I>
NumericMatrix tuples_to_matrix_(List x)
{
  auto p = get_ptr<I>(x);
  NumericMatrix res(p->size(), I);
  transform(begin(*p), end(*p), begin(res), begin(res),
            [&](const array<double, I>& a, double& v)
            {
                auto i = make_strided(&v, p->size());
                copy(begin(a), end(a), i);
                return v;
            });
  return res;
}

// [[Rcpp::export]]
List matrix_to_tuples(const NumericMatrix& x)
{
  switch(x.ncol())
  {
  case 1: return matrix_to_tuples_<1>(x);
  case 2: return matrix_to_tuples_<2>(x);
  case 3: return matrix_to_tuples_<3>(x);
  case 4: return matrix_to_tuples_<4>(x);
  case 5: return matrix_to_tuples_<5>(x);
  case 6: return matrix_to_tuples_<6>(x);
  case 7: return matrix_to_tuples_<7>(x);
  case 8: return matrix_to_tuples_<8>(x);
  case 9: return matrix_to_tuples_<9>(x);
  default: stop("Invalid dimensions");
  }
}

// [[Rcpp::export]]
NumericMatrix tuples_to_matrix(List x)
{
  if (!x.inherits("arrayvec"))
    stop("Expecting arrayvec object");
  switch(as<int>(x["ncol"]))
  {
  case 1: return tuples_to_matrix_<1>(x);
  case 2: return tuples_to_matrix_<2>(x);
  case 3: return tuples_to_matrix_<3>(x);
  case 4: return tuples_to_matrix_<4>(x);
  case 5: return tuples_to_matrix_<5>(x);
  case 6: return tuples_to_matrix_<6>(x);
  case 7: return tuples_to_matrix_<7>(x);
  case 8: return tuples_to_matrix_<8>(x);
  case 9: return tuples_to_matrix_<9>(x);
  default: stop("Invalid dimensions");
  }
}
