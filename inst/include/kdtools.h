#ifndef __MULTIDIM_HPP__
#define __MULTIDIM_HPP__

#include <algorithm>
#include <thread>
#include <tuple>
#include <cmath>

namespace kdtools {

template <size_t I, typename T, typename U>
constexpr U& get(T& x)
{
  return std::get<I>(x);
}

template <typename TupleType>
struct tuple_size
{
  using type = typename std::tuple_size<TupleType>::type;
  using value_type = typename std::tuple_size<TupleType>::value_type;
  static constexpr auto value = std::tuple_size<TupleType>::value;
};

namespace detail {

using std::abs;
using std::next;
using std::size_t;
using std::thread;
using std::distance;
using std::enable_if;
using std::partition;
using std::nth_element;
using std::tuple_element;
using std::iterator_traits;
using std::partition_point;
using std::placeholders::_1;

template <typename T>
T midpos(const T first, const T last)
{
  return next(first, distance(first, last) / 2);
}

template <size_t I>
struct less_nth
{
  template <typename TupleType>
  bool operator()(const TupleType& lhs, const TupleType& rhs)
  {
    return get<I>(lhs) < get<I>(rhs);
  }
};

template<size_t I, typename TupleType>
struct incr_wrap
{
  static constexpr auto value = (I + 1) % tuple_size<TupleType>::value;
};

template<size_t I, typename TupleType>
struct is_not_last
{
  static constexpr auto value = I != tuple_size<TupleType>::value - 1;
};

template<size_t I, typename TupleType>
struct is_last
{
  static constexpr auto value = I == tuple_size<TupleType>::value - 1;
};

template <size_t I, size_t N = 0>
struct kd_less
{
  template <typename TupleType>
  typename enable_if<is_not_last<N, TupleType>::value, bool>::type
  operator()(const TupleType& lhs, const TupleType& rhs) const
  {
    constexpr auto J = incr_wrap<I, TupleType>::value;
    return get<I>(lhs) == get<I>(rhs) ?
      kd_less<J, N + 1>()(lhs, rhs) :
        get<I>(lhs) < get<I>(rhs);
  }
  template <typename TupleType>
  typename enable_if<is_last<N, TupleType>::value, bool>::type
  operator()(const TupleType& lhs, const TupleType& rhs) const
  {
    return get<I>(lhs) < get<I>(rhs);
  }
};

template <typename Pred, size_t I, size_t N = 0>
struct kd_compare
{
  Pred m_pred;
  kd_compare(const Pred& pred) : m_pred(pred) {}
  template <typename TupleType>
  typename enable_if<is_not_last<N, TupleType>::value, bool>::type
  operator()(const TupleType& lhs, const TupleType& rhs) const
  {
    constexpr auto J = incr_wrap<I, TupleType>::value;
    return !m_pred(get<I>(lhs), get<I>(rhs)) &&
      !m_pred(get<I>(rhs), get<I>(lhs)) ?
      kd_compare<Pred, J, N + 1>(m_pred)(lhs, rhs) :
      m_pred(get<I>(lhs), get<I>(rhs));
  }
  template <typename TupleType>
  typename enable_if<is_last<N, TupleType>::value, bool>::type
  operator()(const TupleType& lhs, const TupleType& rhs) const
  {
    return m_pred(get<I>(lhs), get<I>(rhs));
  }
};

template <size_t I, typename Pred, typename T>
T make_kd_compare(const Pred& pred)
{
  return kd_compare<Pred, I>(pred);
}

template <size_t I, typename Iter>
void kd_sort(Iter first, Iter last)
{
  using TupleType = typename iterator_traits<Iter>::value_type;
  constexpr auto J = incr_wrap<I, TupleType>::value;
  if (distance(first, last) > 1)
  {
    auto pred = kd_less<I>();
    auto pivot = midpos(first, last);
    nth_element(first, pivot, last, pred);
    pivot = partition(first, pivot, bind(pred, _1, *pivot));
    kd_sort<J>(next(pivot), last);
    kd_sort<J>(first, pivot);
  }
}

template <size_t I, typename Iter, typename Compare>
void kd_sort(Iter first, Iter last, const Compare& comp)
{
  using TupleType = typename iterator_traits<Iter>::value_type;
  constexpr auto J = (I + 1) % tuple_size<TupleType>::value;
  if (distance(first, last) > 1)
  {
    auto pivot = midpos(first, last);
    auto pred = make_kd_compare<I>(comp);
    nth_element(first,  pivot,  last,  pred);
    pivot = partition(first, pivot, bind(pred, _1, *pivot));
    kd_sort<J>(next(pivot), last, comp);
    kd_sort<J>(first, pivot, comp);
  }
}

template <size_t I, typename Iter>
void kd_sort_threaded(Iter first, Iter last,
                      int max_threads = std::thread::hardware_concurrency(),
                      int thread_depth = 0)
{
  using TupleType = typename iterator_traits<Iter>::value_type;
  constexpr auto J = incr_wrap<I, TupleType>::value;
  if (distance(first, last) > 1)
  {
    auto pred = kd_less<I>();
    auto pivot = midpos(first, last);
    nth_element(first, pivot, last, pred);
    pivot = partition(first, pivot, bind(pred, _1, *pivot));
    if ((1 << thread_depth) < max_threads)
    {
      thread t(kd_sort_threaded<J, Iter>,
               next(pivot), last, max_threads, thread_depth + 1);
      kd_sort_threaded<J>(first, pivot, max_threads, thread_depth + 1);
      t.join();
    }
    else
    {
      kd_sort<J>(next(pivot), last);
      kd_sort<J>(first, pivot);
    }
  }
}

template <size_t I>
struct all_less_
{
  template <typename TupleType>
  typename enable_if<is_not_last<I, TupleType>::value, bool>::type
    operator()(const TupleType& lhs, const TupleType& rhs) const
    {
      return get<I>(lhs) < get<I>(rhs) && all_less_<I + 1>()(lhs, rhs);
    }
  template <typename TupleType>
  typename enable_if<is_last<I, TupleType>::value, bool>::type
    operator()(const TupleType& lhs, const TupleType& rhs) const
    {
      return get<I>(lhs) < get<I>(rhs);
    }
};

template <typename TupleType>
bool all_less(const TupleType& lhs, const TupleType& rhs)
{
  return all_less_<0>()(lhs, rhs);
}

template <size_t I>
struct none_less_
{
  template <typename TupleType>
  typename enable_if<is_not_last<I, TupleType>::value, bool>::type
    operator()(const TupleType& lhs, const TupleType& rhs)
    {
      return get<I>(lhs) >= get<I>(rhs) && none_less_<I + 1>()(lhs, rhs);
    }
  template <typename TupleType>
  typename enable_if<is_last<I, TupleType>::value, bool>::type
    operator()(const TupleType& lhs, const TupleType& rhs)
    {
      return get<I>(lhs) >= get<I>(rhs);
    }
};

template <typename TupleType>
bool none_less(const TupleType& lhs, const TupleType& rhs)
{
  return none_less_<0>()(lhs, rhs);
}

template <size_t I, typename Iter>
Iter find_pivot(Iter first, Iter last)
{
  auto pivot = midpos(first, last);
  return partition_point(first, pivot, bind(less_nth<I>(), _1, *pivot));
}

template <size_t I, typename Iter, typename TupleType>
Iter kd_lower_bound(Iter first, Iter last, const TupleType& value)
{
  constexpr auto J = incr_wrap<I, TupleType>::value;
  if (distance(first, last) > 1)
  {
    auto pivot = find_pivot<I>(first, last);
    if (none_less(*pivot, value))
      return kd_lower_bound<J>(first, pivot, value);
    if (all_less(*pivot, value))
      return kd_lower_bound<J>(next(pivot), last, value);
    auto it = kd_lower_bound<J>(first, pivot, value);
    if (none_less(*it, value)) return it;
    it = kd_lower_bound<J>(next(pivot), last, value);
    if (none_less(*it, value)) return it;
    return last;
  }
  return none_less(*first, value) ? first : last;
}

template <size_t I, typename Iter, typename TupleType>
Iter kd_upper_bound(Iter first, Iter last, const TupleType& value)
{
  constexpr auto J = incr_wrap<I, TupleType>::value;
  if (distance(first, last) > 1)
  {
    auto pivot = find_pivot<I>(first, last);
    if (all_less(value, *pivot))
      return kd_upper_bound<J>(first, pivot, value);
    if (none_less(value, *pivot))
      return kd_upper_bound<J>(next(pivot), last, value);
    auto it = kd_upper_bound<J>(first, pivot, value);
    if (all_less(value, *it)) return it;
    it = kd_upper_bound<J>(next(pivot), last, value);
    if (all_less(value, *it)) return it;
    return last;
  }
  return all_less(value, *first) ? first : last;
}

template <size_t I>
struct sum_of_squares_
{
  template <typename TupleType>
  typename enable_if<is_not_last<I, TupleType>::value, double>::type
    operator()(const TupleType& lhs, const TupleType& rhs) const
    {
      using next_ = sum_of_squares_<I + 1>;
      return std::pow(get<I>(lhs) - get<I>(rhs), 2) + next_()(lhs, rhs);
    }
  template <typename TupleType>
  typename enable_if<is_last<I, TupleType>::value, double>::type
    operator()(const TupleType& lhs, const TupleType& rhs) const
    {
      return std::pow(get<I>(lhs) - get<I>(rhs), 2);
    }
};

template <typename TupleType>
double sum_of_squares(const TupleType& lhs, const TupleType& rhs)
{
  return sum_of_squares_<0>()(lhs, rhs);
}

template <typename TupleType>
double l2dist(const TupleType& lhs, const TupleType& rhs)
{
  return std::sqrt(sum_of_squares(lhs, rhs));
}

template <size_t I, typename Iter, typename TupleType>
Iter kd_nearest_neighbor(Iter first, Iter last, const TupleType& value)
{
  constexpr auto J = incr_wrap<I, TupleType>::value;
  if (distance(first, last) > 1)
  {
    auto pivot = find_pivot<I>(first, last);
    auto search_left = less_nth<I>()(value, *pivot);
    auto search = search_left ?
      kd_nearest_neighbor<J>(first, pivot, value) :
      kd_nearest_neighbor<J>(next(pivot), last, value);
    auto min_dist = l2dist(*pivot, value);
    if (search == last) search = pivot;
    else
    {
      auto sdist = l2dist(*search, value);
      if (sdist < min_dist) min_dist = sdist;
      else search = pivot;
    }
    if (abs(get<I>(value) - get<I>(*pivot)) < min_dist)
    {
      auto s2 = search_left ?
        kd_nearest_neighbor<J>(next(pivot), last, value) :
        kd_nearest_neighbor<J>(first, pivot, value);
      if (s2 != last && l2dist(*s2, value) < min_dist) search = s2;
    }
    return search;
  }
  return first;
}

template <typename TupleType>
bool contains(const TupleType& value,
              const TupleType& lower,
              const TupleType& upper)
{
  return none_less(value, lower) && all_less(value, upper);
}

template <size_t I,
          typename Iter,
          typename TupleType,
          typename OutIter>
void kd_range_query(Iter first, Iter last,
                    const TupleType& lower,
                    const TupleType& upper,
                    OutIter outp)
{
  if (distance(first, last) == 0) return;
  if (distance(first, last) == 1)
  {
    if (contains(*first, lower, upper)) *outp++ = *first;
    return;
  }
  auto pivot = find_pivot<I>(first, last);
  constexpr auto J = incr_wrap<I, TupleType>::value;
  if (contains(*pivot, lower, upper)) *outp++ = *pivot;
  if (!less_nth<I>()(*pivot, lower))
    kd_range_query<J>(first, pivot, lower, upper, outp);
  if (less_nth<I>()(*pivot, upper))
    kd_range_query<J>(next(pivot), last, lower, upper, outp);
  return;
}

}; // namespace detail

template <typename Iter>
void kd_sort(Iter first, Iter last)
{
  detail::kd_sort<0>(first, last);
}

template <typename Iter, typename Compare>
void kd_sort(Iter first, Iter last, const Compare& comp)
{
  detail::kd_sort<0>(first, last, comp);
}

template <typename Iter>
void kd_sort_threaded(Iter first, Iter last)
{
  detail::kd_sort_threaded<0>(first, last);
}

template <typename Iter, typename Value>
Iter kd_lower_bound(Iter first, Iter last, const Value& value)
{
  return detail::kd_lower_bound<0>(first, last, value);
}

template <typename Iter, typename Value>
Iter kd_upper_bound(Iter first, Iter last, const Value& value)
{
  return detail::kd_upper_bound<0>(first, last, value);
}

template <typename Iter, typename TupleType>
bool kd_binary_search(Iter first, Iter last, const TupleType& value)
{
  first = detail::kd_lower_bound<0>(first, last, value);
  return first != last && detail::none_less(value, *first);
}

template <typename Iter, typename Value>
std::pair<Iter, Iter> kd_equal_range(Iter first, Iter last, const Value& value)
{
  return std::make_pair(detail::kd_lower_bound<0>(first, last, value),
                        detail::kd_upper_bound<0>(first, last, value));
}

template <typename Iter, typename TupleType>
Iter kd_nearest_neighbor(Iter first, Iter last, const TupleType& value)
{
  return detail::kd_nearest_neighbor<0>(first, last, value);
}

template <typename Iter,
          typename TupleType,
          typename OutIter>
void kd_range_query(Iter first, Iter last,
                    const TupleType& lower,
                    const TupleType& upper,
                    OutIter outp)
{
  detail::kd_range_query<0>(first, last, lower, upper, outp);
}

template <typename Iter>
void lex_sort(Iter first, Iter last)
{
  std::sort(first, last, detail::kd_less<0>());
}

}; // namespace kdtools

#endif // __MULTIDIM_HPP__
