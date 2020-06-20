#ifndef __TUPLEMAPR_H__
#define __TUPLEMAPR_H__

#include <tuple>
#include <cmath>
#include <array>
#include <string>
#include <utility>
#include <ostream>
#include <type_traits>

#define RUN_TUPLEMAPR_TESTS

namespace keittlab {
namespace tuple {
namespace details {

/*
 * Get type of first element
 */
template<typename... Ts>
using first_of = std::remove_reference_t<
  std::tuple_element_t<0, std::tuple<Ts...>>
>;

/*
 * Index sequence over length of type
 */
template<typename T>
using indices_spanning = std::make_index_sequence<
  std::tuple_size_v<
    std::remove_reference_t<T>
  >
>;

/*
 * Index sequence spanning first tuple type
 */
template<typename... Ts>
using indx_seq_1st_of = indices_spanning<first_of<Ts...>>;

/*
 * Detect void type
 */
template<typename T>
constexpr bool is_void = std::is_same_v<void, T>;

/*
 * Less noisy version of std::forward
 * (c++ needs an unary operator for this)
 */
template<typename T>
constexpr decltype(auto) _(std::remove_reference_t<T>&& t)
{
  return std::forward<T>(t);
}

/*
 * Less noisy version of std::forward
 * (c++ needs an unary operator for this)
 */
template<typename T>
constexpr decltype(auto) _(std::remove_reference_t<T>& t)
{
  return std::forward<T>(t);
}

/*
 * Detect array argument
 */
template<typename>
struct is_std_array : std::false_type {};

/*
* Detect array argument
*/
template<typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};

/*
* Detect array argument
*/
template<typename T>
constexpr bool is_std_array_v = is_std_array<std::decay_t<T>>::value;

/*
* Detect tuple argument
*/
template<typename>
struct is_std_tuple : std::false_type {};

/*
* Detect tuple argument
*/
template<typename... Ts>
struct is_std_tuple<std::tuple<Ts...>> : std::true_type {};

/*
* Detect tuple argument
*/
template<typename T>
constexpr bool is_std_tuple_v = is_std_tuple<std::decay_t<T>>::value;

/*
* Detect pair argument
*/
template<typename>
struct is_std_pair : std::false_type {};

/*
* Detect pair argument
*/
template<typename T, typename U>
struct is_std_pair<std::pair<T, U>> : std::true_type {};

/*
* Detect pair argument
*/
template<typename T>
constexpr bool is_std_pair_v = is_std_pair<std::decay_t<T>>::value;

/**
* Form a tuple of the Ith elements of a set of tuples
*/
template<std::size_t I, typename... Ts>
constexpr decltype(auto) pick(Ts&&... ts) {
 return std::make_tuple(std::get<I>(_<Ts>(ts))...);
}

/*
* Apply a function over the 0th, 1st, 2nd... elements
* of a set of tupples returning a tuple. The pack Ts...
* expands over the set of tuples. The pack Is... expands
* as an integer sequence over the length of the first
* tuple in the set.
*/
template<typename F, std::size_t... Is, typename... Ts>
constexpr decltype(auto) map2tuple_impl(F&& f, std::index_sequence<Is...>, Ts&&... ts) {
 return std::make_tuple(std::apply(_<F>(f), pick<Is>(_<Ts>(ts)...))...);
}

/*
* Applies an invokable and returns an array of results
*/
template<typename F, std::size_t... Is, typename... Ts>
constexpr decltype(auto) map2array_impl(F&& f, std::index_sequence<Is...>, Ts&&... ts) {
 return std::array{std::apply(_<F>(f), pick<Is>(_<Ts>(ts)...))...};
}

/*
* Applies an invokable and returns a pair of results
*/
template<typename F, std::size_t... Is, typename... Ts>
constexpr decltype(auto) map2pair_impl(F&& f, std::index_sequence<Is...>, Ts&&... ts) {
 return std::make_pair(std::apply(_<F>(f), pick<Is>(_<Ts>(ts)...))...);
}

/*
* Applies an invokable and does not return (fold experssions
* are valid targets for pack expansions)
*/
template<typename F, std::size_t... Is, typename... Ts>
constexpr void map2void_impl(F&& f, std::index_sequence<Is...>, Ts&&... ts) {
 (std::apply(_<F>(f), pick<Is>(_<Ts>(ts)...)), ...);
}

/*
* Convenience function to obtain the result of applying only
* to the 0th elements. Used to detect an invokable with a
* void return type.
*/
template<typename F, typename... Ts>
constexpr decltype(auto) map0(F&& f, Ts&&... ts) {
 return std::apply(_<F>(f), pick<0>(_<Ts>(ts)...));
}

} // namespace details

/*
 * Map returning a tuple
 */
template<typename F, typename... Ts>
constexpr decltype(auto) map2tuple(F&& f, Ts&&... ts) {
  using namespace details;
  return map2tuple_impl(_<F>(f), indx_seq_1st_of<Ts...>{}, _<Ts>(ts)...);
}

/*
 * Map returning array
 */
template<typename F, typename... Ts>
constexpr decltype(auto) map2array(F&& f, Ts&&... ts) {
  using namespace details;
  return map2array_impl(_<F>(f), indx_seq_1st_of<Ts...>{}, _<Ts>(ts)...);
}

/*
 * Map returning pair
 */
template<typename F, typename... Ts>
constexpr decltype(auto) map2pair(F&& f, Ts&&... ts) {
  using namespace details;
  return map2pair_impl(_<F>(f), indx_seq_1st_of<Ts...>{}, _<Ts>(ts)...);
}

/*
 * Map with void return
 */
template<typename F, typename... Ts>
constexpr void map2void(F&& f, Ts&&... ts) {
  using namespace details;
  map2void_impl(_<F>(f), indx_seq_1st_of<Ts...>{}, _<Ts>(ts)...);
}

/*
 * Full map function. Automatically handles tuple, array
 * and pair inputs. The return type and the number of elements
 * returned is set by the first type in the Ts... pack. If
 * the invokable F returns void, then void is returned.
 */
template<typename F, typename... Ts>
constexpr decltype(auto) map(F&& f, Ts&&... ts) {
  using namespace details;
  using T = first_of<Ts...>;
  using ret = decltype(map0(_<F>(f), _<Ts>(ts)...));
  if constexpr (is_void<ret>) {
    map2void(_<F>(f), _<Ts>(ts)...);
  } else if constexpr (is_std_pair_v<T>) {
    return map2pair(_<F>(f), _<Ts>(ts)...);
  } else if constexpr (is_std_array_v<T>) {
    return map2array(_<F>(f), _<Ts>(ts)...);
  } else {
    return map2tuple(_<F>(f), _<Ts>(ts)...);
  }
}

namespace details {

/*
 * Since reduction on a single element is
 * a no-op, we reduce the tuple itself if
 * there is only one. Otherwise, we use map
 * to do a parallel reduce over the elements
 * of each tuple.
 */
template<typename F, typename... Ts>
constexpr decltype(auto)
_map(F&& f, Ts&&... ts) {
  using details::_;
  static_assert(sizeof...(ts) > 0);
  if constexpr (sizeof...(ts) > 1) {
    return map(_<F>(f), _<Ts>(ts)...);
  } else {
    return std::apply(_<F>(f), _<Ts>(ts)...);
  }
};

/*
 * Helper avoids static_cast
 */
constexpr double divide(double a, double b) {
  return a / b;
}

} // namespace details

// Unary operations and reductions

/*
 * Not each element
 */
template<typename T>
constexpr decltype(auto)
_not(T&& t) {
  using details::_;
  return map([](auto&& x) {
    return !x;
  }, _<T>(t));
}

/*
 * Sum across tuples or sum
 * of single tuple
 */
template<typename... Ts>
constexpr decltype(auto)
sum(Ts&&... ts) {
  using namespace details;
  return _map([](auto&&... xs) {
    return (xs + ...);
  }, _<Ts>(ts)...);
}

template<typename... Ts>
constexpr decltype(auto)
product(Ts&&... ts) {
  using namespace details;
  return _map([](auto&&... xs) {
    return (xs * ...);
  }, _<Ts>(ts)...);
}

template<typename... Ts>
constexpr decltype(auto)
mean(Ts&&... ts) {
  using namespace details;
  return _map([](auto&&... xs) {
    return divide((xs + ...), sizeof...(xs));
  }, _<Ts>(ts)...);
}

template<typename... Ts>
constexpr decltype(auto)
all_true(Ts&&... ts) {
  using namespace details;
  return _map([](auto&&... xs) {
    return (xs && ...);
  }, _<Ts>(ts)...);
}

template<typename... Ts>
constexpr decltype(auto)
all_false(Ts&&... ts) {
  using namespace details;
  return _map([](auto&&... xs) {
    return (!xs && ...);
  }, _<Ts>(ts)...);
}

template<typename... Ts>
constexpr decltype(auto)
any_true(Ts&&... ts) {
  using namespace details;
  return _map([](auto&&... xs) {
    return (xs || ...);
  }, _<Ts>(ts)...);
}

template<typename... Ts>
constexpr decltype(auto)
any_false(Ts&&... ts) {
  using namespace details;
  return _map([](auto&&... xs) {
    return (!xs || ...);
  }, _<Ts>(ts)...);
}

// binary operations

/*
 * Return constexpr logical tuple indicating
 * where tuple types are not matching
 */
template<typename T, typename U>
constexpr decltype(auto)
is_same(T&& t, U&& u) {
  using details::_;
  return map([](auto&& a, auto&& b) {
    return std::is_same_v<decltype(a), decltype(b)>;
  }, _<T>(t), _<U>(u));
}

template<typename T, typename U>
constexpr decltype(auto)
equal(T&& t, U&& u) {
  using details::_;
  return map([](auto&& a, auto&& b) {
    return a == b;
  }, _<T>(t), _<U>(u));
}

template<typename T, typename U>
constexpr decltype(auto)
all_equal(T&& t, U&& u) {
  using details::_;
  return all_true(equal(_<T>(t), _<U>(u)));
}

template<typename T, typename U>
constexpr decltype(auto)
none_equal(T&& t, U&& u) {
  using details::_;
  return all_false(equal(_<T>(t), _<U>(u)));
}

template<typename T, typename U>
constexpr decltype(auto)
less(T&& t, U&& u) {
  using details::_;
  return map([](auto&& a, auto&& b) {
    return a < b;
  }, _<T>(t), _<U>(u));
}

template<typename T, typename U>
constexpr decltype(auto)
all_less(T&& t, U&& u) {
  using details::_;
  return all_true(less(_<T>(t), _<U>(u)));
}

template<typename T, typename U>
constexpr decltype(auto)
none_less(T&& t, U&& u) {
  using details::_;
  return all_false(less(_<T>(t), _<U>(u)));
}

template<typename T, typename U>
constexpr decltype(auto)
add(T&& t, U&& u) {
  using details::_;
  if constexpr (std::is_arithmetic_v<U>) {
    return map([u](auto&& a) {
      return a + u;
    }, _<T>(t));
  } else {
    return map([](auto&& a, auto&& b) {
      return a + b;
    }, _<T>(t), _<U>(u));
  }
}

template<typename T, typename U>
constexpr decltype(auto)
subtract(T&& t, U&& u) {
  using details::_;
  if constexpr (std::is_arithmetic_v<U>) {
    return map([u](auto&& a) {
      return a - u;
    }, _<T>(t));
  } else {
    return map([](auto&& a, auto&& b) {
      return a - b;
    }, _<T>(t), _<U>(u));
  }
}

template<typename T, typename U>
constexpr decltype(auto)
multiply(T&& t, U&& u) {
  using details::_;
  if constexpr (std::is_arithmetic_v<U>) {
    return map([u](auto&& a) {
      return a * u;
    }, _<T>(t));
  } else {
    return map([](auto&& a, auto&& b) {
      return a * b;
    }, _<T>(t), _<U>(u));
  }
}

template<typename T, typename U>
constexpr decltype(auto)
divide(T&& t, U&& u) {
  using details::_;
  if constexpr (std::is_arithmetic_v<U>) {
    return map([u](auto&& a) {
      return a / u;
    }, _<T>(t));
  } else {
    return map([](auto&& a, auto&& b) {
      return a / b;
    }, _<T>(t), _<U>(u));
  }
}

template<typename T, typename U>
constexpr decltype(auto)
dotprod(T&& t, U&& u) {
  using details::_;
  return sum(multiply(_<T>(t), _<U>(u)));
}

template<typename T, typename U>
constexpr decltype(auto)
hamming(T&& t, U&& u) {
  using details::_;
  return sum(_not(equal(_<T>(t), _<U>(u))));
}

template<typename F, typename T, typename U>
constexpr decltype(auto)
choose(F&& f, T&& t, U&& u) {
  using details::_;
  return map([&f](auto&& a, auto&& b){
    return f() ? a : b;
  }, _<T>(t), _<U>(u));
}

template<typename T, typename U>
constexpr decltype(auto)
wmean(T&& t, U&& u) {
  using details::_;
  return dotprod(_<T>(t), _<U>(u)) / sum(_<U>(u));
}

template<typename T, typename U>
constexpr decltype(auto)
sum_sq_diff(T&& t, U&& u, double exp) {
  using details::_;
  return sum(multiply(subtract(_<T>(t), _<U>(u)), subtract(_<T>(t), _<U>(u))));
}

// not constexpr

template<typename T>
constexpr decltype(auto)
pow(T&& t, double exp) {
  using details::_;
  return map([exp](double base) {
    return std::pow(base, exp);
  }, _<T>(t));
}

template<typename T>
constexpr decltype(auto)
abs(T&& t) {
  using details::_;
  return map([](auto x) {
    return std::abs(x);
  }, _<T>(t));
}

template<typename T>
constexpr decltype(auto)
pnorm(T&& t, double exp) {
  using details::_;
  return std::pow(sum(abs(pow(_<T>(t), exp))), 1 / exp);
}

template<typename T, typename U>
constexpr decltype(auto)
pdist(T&& t, U&& u, double exp) {
  using details::_;
  return pnorm(subtract(_<T>(t), _<U>(u)), exp);
}

template<typename T, typename U>
constexpr decltype(auto)
euclidean_distance(T&& t, U&& u) {
  using details::_;
  return pdist(_<T>(t), _<U>(u), 2);
}

template<typename T, typename U>
constexpr decltype(auto)
manhattan_distance(T&& t, U&& u) {
  using details::_;
  return pdist(_<T>(t), _<U>(u), 1);
}

} // namespace tuple

#ifdef RUN_TUPLEMAPR_TESTS

namespace details {
namespace tests {

constexpr static std::array<double, 3>
a1 = {{1, 2, 3}}, a2 = {{4, 5, 6}}, a3 = {{7, 8, 9}};

constexpr static std::tuple<const char*, double, bool>
  t1{"test", 3.14, false};

constexpr static std::tuple<double, bool, const char*>
  t2{0.33, true, "t2"};

static_assert(tuple::all_false(tuple::is_same(t1, t2)));
static_assert(tuple::all_true(tuple::is_same(a1, a2)));

constexpr static auto t3 = tuple::map([](auto&&... xs){ return true; }, t1, t2, a1);

static_assert(tuple::details::is_std_tuple_v<decltype(t3)>);

constexpr static auto p1 = std::make_pair(1, true);

constexpr static auto p2 = tuple::map([](auto&&...){ return true; }, p1, a1, a2, a3, t1, t2);

static_assert(tuple::details::is_std_pair_v<decltype(p2)>);
static_assert(std::get<0>(p2) && std::get<1>(p2));

constexpr static auto a4 = tuple::map([](auto&&... xs)->double{
  return (0 + ... + xs);
}, a1, a2, a3);

static_assert(tuple::details::is_std_array_v<decltype(a4)>);

static_assert(tuple::sum(a1) == 6);
static_assert(std::get<0>(tuple::sum(a1, a2)) == 5 &&
              std::get<1>(tuple::sum(a1, a2)) == 7 &&
              std::get<2>(tuple::sum(a1, a2)) == 9);

static_assert(tuple::mean(a1) == 2);

constexpr static auto a5 = tuple::map([](double x)->bool{ return x < 4; }, a1);

static_assert(tuple::all_true(a5));
static_assert(tuple::all_false(tuple::_not(a5)));

static_assert(tuple::any_true(a5));
static_assert(!tuple::any_false(a5));

static_assert(tuple::all_less(a1, a2));

static_assert(tuple::sum(tuple::add(a1, a3)) == 30);
static_assert(tuple::sum(tuple::add(a1, 2)) == 12);

static_assert(tuple::sum(tuple::subtract(a1, a3)) == -18);
static_assert(tuple::sum(tuple::subtract(a1, 2)) == 0);

static_assert(tuple::sum(tuple::multiply(a1, a3)) == 50);
static_assert(tuple::sum(tuple::multiply(a1, 2)) == 12);

static_assert(tuple::sum(tuple::divide(a1, a3)) - 0.7261905 < 1e-5);
static_assert(tuple::sum(tuple::divide(a1, 2)) == 3);

static_assert(tuple::dotprod(a1, a3) == 50);

static_assert(tuple::hamming(a1, a2) == 3);

static_assert(tuple::all_equal(tuple::choose([](auto&&...){ return true; }, a1, a2), a1));
static_assert(tuple::all_equal(tuple::choose([](auto&&...){ return false; }, a1, a2), a2));

static_assert(tuple::wmean(a1, a2) - 2.133333 < 1e-5);

} // namespace tests
} // namespace details

#endif // RUN_TUPLEMAPR_TESTS

} // namespace keittlab

#endif // __TUPLEMAPR_H__
