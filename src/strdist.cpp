#include <Rcpp.h>
using namespace Rcpp;

#ifndef NO_CXX17

#include "strdist.h"

#define INIT_TAB {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,                \
                  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,                \
                  2, 0, 0, 0, 0, 0, 0, 0, 0, 0,                \
                  3, 0, 0, 0, 0, 0, 0, 0, 0, 0,                \
                  4, 0, 0, 0, 0, 0, 0, 0, 0, 0,                \
                  5, 0, 0, 0, 0, 0, 0, 0, 0, 0,                \
                  6, 0, 0, 0, 0, 0, 0, 0, 0, 0,                \
                  7, 0, 0, 0, 0, 0, 0, 0, 0, 0,                \
                  8, 0, 0, 0, 0, 0, 0, 0, 0, 0,                \
                  9, 0, 0, 0, 0, 0, 0, 0, 0, 0}

namespace keittlab {
namespace strdist {

int levenshtein(std::string_view s1, std::string_view s2) {
  int n = s1.size() + 1, m = s2.size() + 1;
  thread_local static std::vector<int> tab = INIT_TAB;
  thread_local static int dim = 10;
  if (dim < std::max(n, m)) {
    dim = 2 * std::max(n, m); tab.resize(dim * dim);
    for (int i = 0; i != dim; ++i) tab[i * dim] = i;
    for (int j = 1; j != dim; ++j) tab[j] = j;
  }
  for (int i = 1; i != n; ++i) {
    for (int j = 1; j != m; ++j) {
      int sc = s1[i - 1] == s2[j - 1] ? 0 : 1;
      tab[i * dim + j] = std::min(tab[i * dim + j - 1] + 1,
                                  std::min(tab[(i - 1) * dim + j] + 1,
                                           tab[(i - 1) * dim + j - 1] + sc));
    }
  }
  return tab[(n - 1) * dim + (m - 1)];
}

} // namespace strdist
} // namespace keittlab

#endif // NO_CXX17

//' Compute edit distance
//'
//' @param s1 first character vector
//' @param s2 second character vector
//'
//' @details If neither string is greater than ten characters, there is
//' no allocation. Otherwise, allocation is amortized constant.
//'
//' @return the integer number of edits between strings
//' @export
// [[Rcpp::export]]
int levenshtein(const char* s1, const char* s2) {
#ifdef NO_CXX17
  return NA_INTEGER;
#else
  return keittlab::strdist::levenshtein(s1, s2);
#endif
}
