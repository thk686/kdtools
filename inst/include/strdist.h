#ifndef __STRDIST_H__
#define __STRDIST_H__

#include <string_view>

namespace keittlab {
namespace strdist {

int levenshtein(std::string_view s1, std::string_view s2);

} // namespace strdist
} // namespace keittlab

#endif // __STRDIST_H__
