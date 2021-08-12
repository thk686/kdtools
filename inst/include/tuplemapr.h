#ifndef TUPLEMAPR_H
#define TUPLEMAPR_H

#include <tuple>
#include <cmath>
#include <array>
#include <string>
#include <utility>
#include <ostream>
#include <type_traits>

// #define RUN_TUPLEMAPR_TESTS

namespace keittlab {
namespace tuple {
namespace detail {

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
 return std::make_tuple(std::get<I>(std::forward<Ts>(ts))...);
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
 return std::make_tuple(std::apply(std::forward<F>(f), pick<Is>(std::forward<Ts>(ts)...))...);
}

/*
* Applies an invokable and returns an array of results
*/
template<typename F, std::size_t... Is, typename... Ts>
constexpr decltype(auto) map2array_impl(F&& f, std::index_sequence<Is...>, Ts&&... ts) {
 return std::array{std::apply(std::forward<F>(f), pick<Is>(std::forward<Ts>(ts)...))...};
}

/*
* Applies an invokable and returns a pair of results
*/
template<typename F, std::size_t... Is, typename... Ts>
constexpr decltype(auto) map2pair_impl(F&& f, std::index_sequence<Is...>, Ts&&... ts) {
 return std::make_pair(std::apply(std::forward<F>(f), pick<Is>(std::forward<Ts>(ts)...))...);
}

/*
* Applies an invokable and does not return (fold experssions
* are valid targets for pack expansions)
*/
template<typename F, std::size_t... Is, typename... Ts>
constexpr void map2void_impl(F&& f, std::index_sequence<Is...>, Ts&&... ts) {
 (std::apply(std::forward<F>(f), pick<Is>(std::forward<Ts>(ts)...)), ...);
}

/*
* Convenience function to obtain the result of applying only
* to the 0th elements. Used to detect an invokable with a
* void return type.
*/
template<typename F, typename... Ts>
constexpr decltype(auto) map0(F&& f, Ts&&... ts) {
 return std::apply(std::forward<F>(f), pick<0>(std::forward<Ts>(ts)...));
}

/*
 * Helper less busy than writing static_cast...
 */
constexpr double divide(double a, double b) {
  return a / b;
}

} // namespace detail

/*
 * Map returning a tuple
 */
template<typename F, typename... Ts>
constexpr decltype(auto) map2tuple(F&& f, Ts&&... ts) {
  return map2tuple_impl(std::forward<F>(f), detail::indx_seq_1st_of<Ts...>{}, std::forward<Ts>(ts)...);
}

/*
 * Map returning array
 */
template<typename F, typename... Ts>
constexpr decltype(auto) map2array(F&& f, Ts&&... ts) {
  return map2array_impl(std::forward<F>(f), detail::indx_seq_1st_of<Ts...>{}, std::forward<Ts>(ts)...);
}

/*
 * Map returning pair
 */
template<typename F, typename... Ts>
constexpr decltype(auto) map2pair(F&& f, Ts&&... ts) {
  return map2pair_impl(std::forward<F>(f), detail::indx_seq_1st_of<Ts...>{}, std::forward<Ts>(ts)...);
}

/*
 * Map with void return
 */
template<typename F, typename... Ts>
constexpr void map2void(F&& f, Ts&&... ts) {
  map2void_impl(std::forward<F>(f), detail::indx_seq_1st_of<Ts...>{}, std::forward<Ts>(ts)...);
}

/*
 * Full map function. Automatically handles tuple, array
 * and pair inputs. The return type and the number of elements
 * returned is set by the first type in the Ts... pack. If
 * the invokable F returns void, then void is returned.
 */
template<typename F, typename... Ts>
constexpr decltype(auto) map(F&& f, Ts&&... ts) {
  using namespace detail;
  using T = first_of<Ts...>;
  using ret = decltype(map0(std::forward<F>(f), std::forward<Ts>(ts)...));
  if constexpr (is_void<ret>) {
    map2void(std::forward<F>(f), std::forward<Ts>(ts)...);
  } else if constexpr (is_std_pair_v<T>) {
    return map2pair(std::forward<F>(f), std::forward<Ts>(ts)...);
  } else if constexpr (is_std_array_v<T>) {
    return map2array(std::forward<F>(f), std::forward<Ts>(ts)...);
  } else {
    return map2tuple(std::forward<F>(f), std::forward<Ts>(ts)...);
  }
}

/*
 * Since reduction on a single element is
 * a no-op, we reduce the tuple itself if
 * there is only one. Otherwise, we use map
 * to do a parallel reduce over the elements
 * of each tuple.
 */
template<typename F, typename... Ts>
constexpr decltype(auto)
map_reduce(F&& f, Ts&&... ts) {
  static_assert(sizeof...(ts) > 0);
  if constexpr (sizeof...(ts) > 1) {
    return map(std::forward<F>(f), std::forward<Ts>(ts)...);
  } else {
    return std::apply(std::forward<F>(f), std::forward<Ts>(ts)...);
  }
}

// Unary operations and reductions

/*
 * Not each element
 */
template<typename T>
constexpr decltype(auto)
_not(T&& t) {
  return map([](auto&& x) {
    return !x;
  }, std::forward<T>(t));
}

/*
 * Sum across tuples or sum
 * of single tuple
 */
template<typename... Ts>
constexpr decltype(auto)
sum(Ts&&... ts) {
  return map_reduce([](auto&&... xs) {
    return (xs + ...);
  }, std::forward<Ts>(ts)...);
}

template<typename... Ts>
constexpr decltype(auto)
product(Ts&&... ts) {
  return map_reduce([](auto&&... xs) {
    return (xs * ...);
  }, std::forward<Ts>(ts)...);
}

template<typename... Ts>
constexpr decltype(auto)
mean(Ts&&... ts) {
  return map_reduce([](auto&&... xs) {
    return divide((xs + ...), sizeof...(xs));
  }, std::forward<Ts>(ts)...);
}

template<typename... Ts>
constexpr decltype(auto)
all_true(Ts&&... ts) {
  return map_reduce([](auto&&... xs) {
    return (xs && ...);
  }, std::forward<Ts>(ts)...);
}

template<typename... Ts>
constexpr decltype(auto)
all_false(Ts&&... ts) {
  return all_true(_not(std::forward<Ts>(ts)...));
}

template<typename... Ts>
constexpr decltype(auto)
any_true(Ts&&... ts) {
  return map_reduce([](auto&&... xs) {
    return (xs || ...);
  }, std::forward<Ts>(ts)...);
}

template<typename... Ts>
constexpr decltype(auto)
any_false(Ts&&... ts) {
  return any_true(_not(std::forward<Ts>(ts)...));
}

// binary operations

/*
 * Return constexpr logical tuple indicating
 * where tuple types are not matching
 */
template<typename T, typename U>
constexpr decltype(auto)
is_same(T&& t, U&& u) {
  return map([](auto&& a, auto&& b) {
    return std::is_same_v<decltype(a), decltype(b)>;
  }, std::forward<T>(t), std::forward<U>(u));
}

template<typename T, typename U>
constexpr decltype(auto)
equal(T&& t, U&& u) {
  return map([](auto&& a, auto&& b) {
    return a == b;
  }, std::forward<T>(t), std::forward<U>(u));
}

template<typename T, typename U>
constexpr decltype(auto)
all_equal(T&& t, U&& u) {
  return all_true(equal(std::forward<T>(t), std::forward<U>(u)));
}

template<typename T, typename U>
constexpr decltype(auto)
none_equal(T&& t, U&& u) {
  return all_false(equal(std::forward<T>(t), std::forward<U>(u)));
}

template<typename T, typename U>
constexpr decltype(auto)
less(T&& t, U&& u) {
  return map([](auto&& a, auto&& b) {
    return a < b;
  }, std::forward<T>(t), std::forward<U>(u));
}

template<typename T, typename U>
constexpr decltype(auto)
all_less(T&& t, U&& u) {
  return all_true(less(std::forward<T>(t), std::forward<U>(u)));
}

template<typename T, typename U>
constexpr decltype(auto)
none_less(T&& t, U&& u) {
  return all_false(less(std::forward<T>(t), std::forward<U>(u)));
}

template<typename T, typename U>
constexpr decltype(auto)
add(T&& t, U&& u) {
  if constexpr (std::is_arithmetic_v<U>) {
    return map([u](auto&& a) {
      return a + u;
    }, std::forward<T>(t));
  } else {
    return map([](auto&& a, auto&& b) {
      return a + b;
    }, std::forward<T>(t), std::forward<U>(u));
  }
}

template<typename T, typename U>
constexpr decltype(auto)
subtract(T&& t, U&& u) {
  if constexpr (std::is_arithmetic_v<U>) {
    return map([u](auto&& a) {
      return a - u;
    }, std::forward<T>(t));
  } else {
    return map([](auto&& a, auto&& b) {
      return a - b;
    }, std::forward<T>(t), std::forward<U>(u));
  }
}

template<typename T, typename U>
constexpr decltype(auto)
multiply(T&& t, U&& u) {
  if constexpr (std::is_arithmetic_v<U>) {
    return map([u](auto&& a) {
      return a * u;
    }, std::forward<T>(t));
  } else {
    return map([](auto&& a, auto&& b) {
      return a * b;
    }, std::forward<T>(t), std::forward<U>(u));
  }
}

template<typename T, typename U>
constexpr decltype(auto)
divide(T&& t, U&& u) {
  if constexpr (std::is_arithmetic_v<U>) {
    return map([u](auto&& a) {
      return a / u;
    }, std::forward<T>(t));
  } else {
    return map([](auto&& a, auto&& b) {
      return a / b;
    }, std::forward<T>(t), std::forward<U>(u));
  }
}

template<typename T, typename U>
constexpr decltype(auto)
dotprod(T&& t, U&& u) {
  return sum(multiply(std::forward<T>(t), std::forward<U>(u)));
}

template<typename T, typename U>
constexpr decltype(auto)
hamming(T&& t, U&& u) {
  return sum(_not(equal(std::forward<T>(t), std::forward<U>(u))));
}

template<typename F, typename T, typename U>
constexpr decltype(auto)
choose(F&& f, T&& t, U&& u) {
  return map([&f](auto&& a, auto&& b){
    return f() ? a : b;
  }, std::forward<T>(t), std::forward<U>(u));
}

template<typename T, typename U>
constexpr decltype(auto)
wmean(T&& t, U&& u) {
  return dotprod(std::forward<T>(t), std::forward<U>(u)) / sum(std::forward<U>(u));
}

template<typename T, typename U>
constexpr decltype(auto)
sum_sq_diff(T&& t, U&& u, double exp) {
  return sum(multiply(subtract(std::forward<T>(t), std::forward<U>(u)), subtract(std::forward<T>(t), std::forward<U>(u))));
}

// not constexpr

template<typename T>
constexpr decltype(auto)
pow(T&& t, double exp) {
  return map([exp](double base) {
    return std::pow(base, exp);
  }, std::forward<T>(t));
}

template<typename T>
constexpr decltype(auto)
abs(T&& t) {
  return map([](auto x) {
    return std::abs(x);
  }, std::forward<T>(t));
}

template<typename T>
constexpr decltype(auto)
pnorm(T&& t, double exp) {
  return std::pow(sum(abs(pow(std::forward<T>(t), exp))), 1 / exp);
}

template<typename T, typename U>
constexpr decltype(auto)
pdist(T&& t, U&& u, double exp) {
  return pnorm(subtract(std::forward<T>(t), std::forward<U>(u)), exp);
}

template<typename T, typename U>
constexpr decltype(auto)
euclidean_distance(T&& t, U&& u) {
  return pdist(std::forward<T>(t), std::forward<U>(u), 2);
}

template<typename T, typename U>
constexpr decltype(auto)
manhattan_distance(T&& t, U&& u) {
  return pdist(std::forward<T>(t), std::forward<U>(u), 1);
}

} // namespace tuple

#ifdef RUN_TUPLEMAPR_TESTS

namespace detail {
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

static_assert(tuple::detail::is_std_tuple_v<decltype(t3)>);

constexpr static auto p1 = std::make_pair(1, true);

constexpr static auto p2 = tuple::map([](auto&&...){ return true; }, p1, a1, a2, a3, t1, t2);

static_assert(tuple::detail::is_std_pair_v<decltype(p2)>);
static_assert(std::get<0>(p2) && std::get<1>(p2));

constexpr static auto a4 = tuple::map([](auto&&... xs)->double{
  return (0 + ... + xs);
}, a1, a2, a3);

static_assert(tuple::detail::is_std_array_v<decltype(a4)>);

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
} // namespace detail

#endif // RUN_TUPLEMAPR_TESTS

} // namespace keittlab

#endif // TUPLEMAPR_H
