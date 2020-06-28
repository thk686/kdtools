#include <Rcpp.h>
using namespace Rcpp;

#include "strdist.h"

namespace keittlab {
namespace strdist {

int levenshtein(std::string_view s1, std::string_view s2) {
  int n = s1.size() + 1, m = s2.size() + 1;
  static std::vector<int> tab;
  static int lim = 0;
  if (lim < std::max(n, m)) {
    lim = 2 * std::max(n, m); tab.resize(lim * lim);
    for (int i = 0; i != lim; ++i) tab[i * lim] = i;
    for (int j = 1; j != lim; ++j) tab[j] = j;
  }
  for (int i = 1; i != n; ++i) {
    for (int j = 1; j != m; ++j) {
      int sc = s1[i - 1] == s2[j - 1] ? 0 : 1;
      tab[i * lim + j] = std::min(tab[i * lim + j - 1] + 1,
                                  std::min(tab[(i - 1) * lim + j] + 1,
                                           tab[(i - 1) * lim + j - 1] + sc));
    }
  }
  return tab[(n - 1) * lim + (m - 1)];
}

}; // namespace strdist
}; // namespace keittlab

// [[Rcpp::export]]
int levenshtein(const char* s1, const char* s2) {
  return keittlab::strdist::levenshtein(s1, s2);
}
