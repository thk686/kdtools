#ifndef __STRDIST_H__
#define __STRDIST_H__

#ifndef NO_CXX17

#include <string_view>

// Export for other packages to use

namespace keittlab {
namespace strdist {

int levenshtein(std::string_view s1, std::string_view s2);

} // namespace strdist
} // namespace keittlab

#endif // NO_CXX17

#endif // __STRDIST_H__
