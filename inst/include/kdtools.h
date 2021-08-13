// Copyright Timothy H. Keitt 2020

#ifndef KDTOOLS_H
#define KDTOOLS_H

#ifndef NO_CXX17 // CRAN

#include <functional>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <utility>
#include <thread>
#include <vector>
#include <limits>
#include <queue>
#include <tuple>
#include <cmath>

// #define USE_CIRCULAR_LEXICOGRAPHIC_COMPARE
// #define NO_TUPLEMAPR

#ifndef NO_TUPLEMAPR
#include "tuplemapr.h"
#endif // NO_TUPLEMAPR

namespace keittlab {
namespace kdtools {

template <typename T> static constexpr auto
  ndim = std::tuple_size_v<typename std::remove_pointer_t<T>>;

template <typename T>
double scalar_dist(const T& lhs, const T& rhs)
{
  return std::abs(lhs - rhs);
}

namespace detail {

using std::abs;
using std::get;
using std::next;
using std::prev;
using std::pair;
using std::sort;
using std::size_t;
using std::thread;
using std::vector;
using std::is_same;
using std::distance;
using std::pop_heap;
using std::push_heap;
using std::make_heap;
using std::enable_if;
using std::transform;
using std::partition;
using std::nth_element;
using std::is_pointer_v;
using std::tuple_size_v;
using std::numeric_limits;
using std::priority_queue;
using std::is_partitioned;
using std::remove_pointer;
using std::iterator_traits;
using std::partition_point;
using std::placeholders::_1;

template <typename T>
T middle_of(const T first, const T last)
{
  return next(first, distance(first, last) / 2);
}

template <size_t I>
struct less_nth
{
  template <typename T>
  bool operator()(const T& lhs, const T& rhs)
  {
    if constexpr (is_pointer_v<T>)
      return get<I>(*lhs) < get<I>(*rhs);
    else
      return get<I>(lhs) < get<I>(rhs);
  }
  template <typename T, typename U>
  bool operator()(const pair<T, U>& lhs, const pair<T, U>& rhs)
  {
    if constexpr (is_pointer_v<T>)
      return get<I>(*lhs.first) < get<I>(*rhs.first);
    else
      return get<I>(lhs.first) < get<I>(rhs.first);
  }
};

template <size_t I>
struct less_radius_nth
{
  template <typename T, typename U>
  bool operator()(const T& lhs, const T& rhs, const U radius)
  {
    if constexpr (is_pointer_v<T>)
      return scalar_dist(get<I>(*lhs), get<I>(*rhs)) < radius;
    else
      return scalar_dist(get<I>(lhs), get<I>(rhs)) < radius;
  }
  template <typename T, typename U, typename V>
  bool operator()(const pair<T, U>& lhs, const pair<T, U>& rhs, const V radius)
  {
    if constexpr (is_pointer_v<T>)
      return scalar_dist(get<I>(*lhs.first), get<I>(*rhs.first)) < radius;
    else
      return scalar_dist(get<I>(lhs.first), get<I>(rhs.first)) < radius;
  }
};

template <size_t I>
struct equal_nth
{
  template <typename T>
  bool operator()(const T& lhs, const T& rhs)
  {
    if constexpr (is_pointer_v<T>)
      return get<I>(*lhs) == get<I>(*rhs);
    else
      return get<I>(lhs) == get<I>(rhs);
  }
  template <typename T, typename U>
  bool operator()(const pair<T, U>& lhs, const pair<T, U>& rhs)
  {
    if constexpr (is_pointer_v<T>)
      return get<I>(*lhs.first) == get<I>(*rhs.first);
    else
      return get<I>(lhs.first) == get<I>(rhs.first);
  }
};

template <typename Pred, size_t I>
struct pred_nth
{
  Pred m_pred;
  pred_nth(const Pred& pred) : m_pred(pred) {}
  template <typename T>
  bool operator()(const T& lhs, const T& rhs)
  {
    if constexpr (is_pointer_v<T>)
      return m_pred(get<I>(*lhs), get<I>(*rhs));
    else
      return m_pred(get<I>(lhs), get<I>(rhs));
  }
  template <typename T, typename U>
  bool operator()(const pair<T, U>& lhs, const pair<T, U>& rhs)
  {
    if constexpr (is_pointer_v<T>)
      return m_pred(get<I>(*lhs.first), get<I>(*rhs.first));
    else
      return m_pred(get<I>(lhs.first), get<I>(rhs.first));
  }
};

template <size_t I, typename Pred>
pred_nth<Pred, I> make_pred_nth(const Pred& pred)
{
  return pred_nth<Pred, I>(pred);
}

template <size_t I, typename T>
double dist_nth(const T& lhs, const T& rhs)
{
  if constexpr (is_pointer_v<T>)
    return scalar_dist(get<I>(*lhs), get<I>(*rhs));
  else
    return scalar_dist(get<I>(lhs), get<I>(rhs));
}

template <size_t I, typename T, typename U>
double dist_nth(const pair<T, U>& lhs, const pair<T, U>& rhs)
{
  if constexpr (is_pointer_v<T>)
    return scalar_dist(get<I>(*lhs.first), get<I>(*rhs.first));
  else
    return scalar_dist(get<I>(lhs.first), get<I>(rhs.first));
}

template <size_t I, typename T>
static constexpr auto next_dim = (I + 1) % ndim<T>;

template<size_t I, typename T>
static constexpr auto is_last = (I == (ndim<T> - 1));

#ifdef USE_CIRCULAR_LEXICOGRAPHIC_COMPARE

template <size_t I, size_t K = 0>
struct kd_less
{
  template <typename T>
  bool operator()(const T& lhs, const T& rhs) const
  {
    if constexpr (is_last<K, T>) {
      return less_nth<I>()(lhs, rhs);
    } else {
      constexpr auto J = next_dim<I, T>;
      return equal_nth<I>()(lhs, rhs) ?
        kd_less<J, K + 1>()(lhs, rhs) :
          less_nth<I>()(lhs, rhs);
    }
  }
};

template <typename Pred, size_t I, size_t K = 0>
struct kd_compare
{
  Pred m_pred;
  kd_compare(const Pred& pred) : m_pred(pred) {}
  template <typename T>
  bool operator()(const T& lhs, const T& rhs) const
  {
    if constexpr (is_last<K, T>) {
      return make_pred_nth<I>(m_pred)(lhs, rhs);
    } else {
      constexpr auto J = next_dim<I, T>;
      auto pred = make_pred_nth<I>(m_pred);
      return !pred(lhs, rhs) && !pred(rhs, lhs) ?
        kd_compare<Pred, J, K + 1>(m_pred)(lhs, rhs) : pred(lhs, rhs);
    }
  }
};

#else // USE_CIRCULAR_LEXICOGRAPHIC_COMPARE

template <size_t I>
struct kd_less
{
  template <typename T>
  bool operator()(const T& lhs, const T& rhs) const {
    return less_nth<I>()(lhs, rhs);
  }
};

template <typename Pred, size_t I, size_t K = 0>
struct kd_compare
{
  Pred m_pred;
  kd_compare(const Pred& pred) : m_pred(pred) {}
  template <typename T>
  bool operator()(const T& lhs, const T& rhs) const
  {
      return make_pred_nth<I>(m_pred)(lhs, rhs);
  }
};

#endif // USE_CIRCULAR_LEXICOGRAPHIC_COMPARE

template <size_t I, typename Pred>
kd_compare<Pred, I> make_kd_compare(const Pred& pred)
{
  return kd_compare<Pred, I>(pred);
}

template <typename T>
using iter_value_t = typename iterator_traits<T>::value_type;

template <typename Iter, typename Pred>
Iter median_part(Iter first, Iter last, Pred pred)
{
  auto pivot = middle_of(first, last);
  nth_element(first, pivot, last, pred);
  return pivot;
}

template <size_t I, typename Iter>
void kd_sort(Iter first, Iter last)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>;
  if (distance(first, last) > 1)
  {
    auto pred = kd_less<I>();
    auto pivot = median_part(first, last, pred);
    kd_sort<J>(next(pivot), last);
    kd_sort<J>(first, pivot);
  }
}

template <size_t I, typename Iter, typename Compare>
void kd_sort(Iter first, Iter last, const Compare& comp)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>;
  if (distance(first, last) > 1)
  {
    auto pred = make_kd_compare<I>(comp);
    auto pivot = median_part(first, last, pred);
    kd_sort<J>(next(pivot), last, comp);
    kd_sort<J>(first, pivot, comp);
  }
}

template <size_t I, typename Iter>
void kd_sort_threaded(Iter first, Iter last,
                      int max_threads = std::thread::hardware_concurrency(),
                      int thread_depth = 1)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>;
  if (distance(first, last) > 1)
  {
    auto pred = kd_less<I>();
    auto pivot = median_part(first, last, pred);
    if ((1 << thread_depth) <= max_threads)
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

template <size_t I, typename Iter, typename Compare>
void kd_sort_threaded(Iter first, Iter last, const Compare& comp,
                      int max_threads = std::thread::hardware_concurrency(),
                      int thread_depth = 1)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>;
  if (distance(first, last) > 1)
  {
    auto pred = make_kd_compare<I>(comp);
    auto pivot = median_part(first, last, pred);
    if ((1 << thread_depth) <= max_threads)
    {
      thread t(kd_sort_threaded<J, Iter, Compare>,
               next(pivot), last, comp, max_threads, thread_depth + 1);
      kd_sort_threaded<J>(first, pivot, comp, max_threads, thread_depth + 1);
      t.join();
    }
    else
    {
      kd_sort<J>(next(pivot), last, comp);
      kd_sort<J>(first, pivot, comp);
    }
  }
}

template <typename Iter, typename Pred>
bool check_partition(Iter first, Iter pivot, Iter last, Pred pred)
{
  while (first != pivot) if (pred(*pivot, *first++)) return false;
  while (first != last) if (pred(*first++, *pivot)) return false;
  return true;
}

template <size_t I, typename Iter>
bool kd_is_sorted(Iter first, Iter last)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>;
  if (distance(first, last) < 2) return true;
  auto pred = kd_less<I>();
  auto pivot = middle_of(first, last);
  return check_partition(first, pivot, last, pred) &&
    kd_is_sorted<J>(first, pivot) &&
    kd_is_sorted<J>(next(pivot), last);
}

template <size_t I, typename Iter, typename Compare>
bool kd_is_sorted(Iter first, Iter last, const Compare& comp)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>;
  if (distance(first, last) < 2) return true;
  auto pred = make_kd_compare<I>(comp);
  auto pivot = middle_of(first, last);
  return check_partition(first, pivot, last, pred) &&
    kd_is_sorted<J>(first, pivot, comp) &&
    kd_is_sorted<J>(next(pivot), last, comp);
}

template <size_t I, typename Iter>
bool kd_is_sorted_threaded(Iter first, Iter last,
                           int max_threads = std::thread::hardware_concurrency(),
                           int thread_depth = 1)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>;
  if (distance(first, last) < 2) return true;
  auto pred = kd_less<I>();
  auto pivot = middle_of(first, last);
  if (check_partition(first, pivot, last, pred)) {
    if ((1 << thread_depth) <= max_threads)
    {
      bool res_left, res_right;
      thread t([=, &res_left](){
        res_left = kd_is_sorted_threaded<J>(first, pivot, max_threads, thread_depth + 1);
      });
      res_right = kd_is_sorted_threaded<J>(next(pivot), last, max_threads, thread_depth + 1);
      t.join();
      return res_left && res_right;
    }
    else
    {
      return kd_is_sorted<J>(first, pivot) && kd_is_sorted<J>(next(pivot), last);
    }
  } else {
    return false;
  }
}

template <size_t I, typename Iter, typename Compare>
bool kd_is_sorted_threaded(Iter first, Iter last, const Compare& comp,
                           int max_threads = std::thread::hardware_concurrency(),
                           int thread_depth = 1)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>;
  if (distance(first, last) < 2) return true;
  auto pred = make_kd_compare<I>(comp);
  auto pivot = middle_of(first, last);
  if (check_partition(first, pivot, last, pred)) {
    if ((1 << thread_depth) <= max_threads)
    {
      bool res_left, res_right;
      thread t([=, &res_left](){
        res_left = kd_is_sorted_threaded<J>(first, pivot, comp, max_threads, thread_depth + 1);
      });
      res_right = kd_is_sorted_threaded<J>(next(pivot), last, comp, max_threads, thread_depth + 1);
      t.join();
      return res_left && res_right;
    }
    else
    {
      return kd_is_sorted<J>(first, pivot, comp) && kd_is_sorted<J>(next(pivot), last, comp);
    }
  } else {
    return false;
  }
}

#ifdef NO_TUPLEMAPR

// leftover from c++11 version
// tuplemapr can fail on some older libc++ versions

template <size_t I>
struct all_less_
{
  template <typename T>
  bool operator()(const T& lhs, const T& rhs) const
  {
    if constexpr (is_last<I, T>) {
      return less_nth<I>()(lhs, rhs);
    } else {
      return less_nth<I>()(lhs, rhs) && all_less_<I + 1>()(lhs, rhs);
    }
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
  template <typename T>
  bool operator()(const T& lhs, const T& rhs)
  {
    if constexpr (is_last<I, T>) {
      return !less_nth<I>()(lhs, rhs);
    } else {
      return !less_nth<I>()(lhs, rhs) && none_less_<I + 1>()(lhs, rhs);
    }
  }
};

template <typename TupleType>
bool none_less(const TupleType& lhs, const TupleType& rhs)
{
  return none_less_<0>()(lhs, rhs);
}

template <size_t I>
struct sum_of_squares_
{
  template <typename TupleType>
  double operator()(const TupleType& lhs, const TupleType& rhs) const
  {
    if constexpr (is_last<I, TupleType>) {
      return std::pow(dist_nth<I>(rhs, lhs), 2);
    } else {
      using next_ = sum_of_squares_<I + 1>;
      return std::pow(dist_nth<I>(rhs, lhs), 2) + next_()(lhs, rhs);
    }
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

template <size_t I>
struct p_sum_
{
  template <typename TupleType>
  double operator()(const TupleType& lhs,
                  const TupleType& rhs,
                  double p) const
  {
    if constexpr (is_last<I, TupleType>) {
      return std::pow(dist_nth<I>(rhs, lhs), p);
    } else {
      using next_ = p_sum_<I + 1>;
      return std::pow(dist_nth<I>(rhs, lhs), p) + next_()(lhs, rhs, p);
    }
  }
};

template <typename TupleType>
double p_sum(const TupleType& lhs, const TupleType& rhs, double p)
{
  return p_sum_<0>()(lhs, rhs, p);
}

template <typename TupleType>
double pdist(const TupleType& lhs, const TupleType& rhs, double p)
{
  return std::pow(p_sum(lhs, rhs, p), 1 / p);
}

#else // NO_TUPLEMAPR

using tuple::all_less;
using tuple::none_less;
using tuple::pdist;
using l2dist = tuple::euclidean_distance;

#endif // NO_TUPLEMAPR

template <size_t I, typename Iter, typename TupleType>
Iter kd_lower_bound(Iter first, Iter last, const TupleType& key)
{
  constexpr auto J = next_dim<I, TupleType>;
  if (distance(first, last) > 1)
  {
    auto pivot = middle_of(first, last);
    if (none_less(*pivot, key))
      return kd_lower_bound<J>(first, pivot, key);
    if (all_less(*pivot, key))
      return kd_lower_bound<J>(next(pivot), last, key);
    auto it = kd_lower_bound<J>(first, pivot, key);
    if (it != last && none_less(*it, key)) return it;
    it = kd_lower_bound<J>(next(pivot), last, key);
    if (it != last && none_less(*it, key)) return it;
    return last;
  }
  return first != last && none_less(*first, key) ? first : last;
}

template <size_t I, typename Iter, typename TupleType>
Iter kd_upper_bound(Iter first, Iter last, const TupleType& key)
{
  constexpr auto J = next_dim<I, TupleType>;
  if (distance(first, last) > 1)
  {
    auto pivot = middle_of(first, last);
    if (all_less(key, *pivot))
      return kd_upper_bound<J>(first, pivot, key);
    if (none_less(key, *pivot))
      return kd_upper_bound<J>(next(pivot), last, key);
    auto it = kd_upper_bound<J>(first, pivot, key);
    if (it != last && all_less(key, *it)) return it;
    it = kd_upper_bound<J>(next(pivot), last, key);
    if (it != last && all_less(key, *it)) return it;
    return last;
  }
  return first != last && all_less(key, *first) ? first : last;
}

template <size_t I, typename Iter, typename TupleType>
Iter kd_nearest_neighbor(Iter first, Iter last, const TupleType& key)
{
  constexpr auto J = next_dim<I, TupleType>;
  if (distance(first, last) > 1)
  {
    auto pivot = middle_of(first, last);
    if (equal_nth<I>()(*pivot, key)) {
      auto left_res = kd_nearest_neighbor<J>(first, pivot, key);
      auto right_res = kd_nearest_neighbor<J>(next(pivot), last, key);
      if (l2dist(*right_res, key) < l2dist(*left_res, key)) {
        return right_res;
      } else {
        return left_res;
      }
    } else {
      auto search_left = less_nth<I>()(key, *pivot);
      auto search = search_left ?
        kd_nearest_neighbor<J>(first, pivot, key) :
          kd_nearest_neighbor<J>(next(pivot), last, key);
      auto min_dist = l2dist(*pivot, key);
      if (search == last) {
        search = pivot;
      } else {
        auto sdist = l2dist(*search, key);
        if (sdist < min_dist) min_dist = sdist;
        else search = pivot;
      }
      if (dist_nth<I>(key, *pivot) < min_dist) {
      auto s2 = search_left ?
        kd_nearest_neighbor<J>(next(pivot), last, key) :
          kd_nearest_neighbor<J>(first, pivot, key);
      if (s2 != last && l2dist(*s2, key) < min_dist) search = s2;
      }
     return search;
    }
  }
  return first;
}

template <size_t I, typename Iter, typename TupleType>
Iter kd_nearest_neighbor(Iter first, Iter last,
                         const TupleType& key,
                         double p)
{
  constexpr auto J = next_dim<I, TupleType>;
  if (distance(first, last) > 1)
  {
    auto pivot = middle_of(first, last);
    if (equal_nth<I>()(*pivot, key)) {
      auto left_res = kd_nearest_neighbor<J>(first, pivot, key);
      auto right_res = kd_nearest_neighbor<J>(next(pivot), last, key);
      if (pdist(*right_res, key, p) < pdist(*left_res, key, p)) {
        return right_res;
      } else {
        return left_res;
      }
    } else {
      auto search_left = less_nth<I>()(key, *pivot);
      auto search = search_left ?
      kd_nearest_neighbor<J>(first, pivot, key) :
        kd_nearest_neighbor<J>(next(pivot), last, key);
      auto min_dist = pdist(*pivot, key, p);
      if (search == last) {
        search = pivot;
      } else {
        auto sdist = pdist(*search, key, p);
        if (sdist < min_dist) min_dist = sdist;
        else search = pivot;
      }
      if (dist_nth<I>(key, *pivot) < min_dist) {
        auto s2 = search_left ?
        kd_nearest_neighbor<J>(next(pivot), last, key) :
        kd_nearest_neighbor<J>(first, pivot, key);
        if (s2 != last && pdist(*s2, key, p) < min_dist) search = s2;
      }
      return search;
    }
  }
  return first;
}

template <typename TupleType>
bool within(const TupleType& key,
            const TupleType& lower,
            const TupleType& upper)
{
  return none_less(key, lower) && all_less(key, upper);
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
  if (distance(first, last) > 32) {
    auto pred = less_nth<I>();
    auto pivot = middle_of(first, last);
    constexpr auto J = next_dim<I, TupleType>;
    if (within(*pivot, lower, upper)) *outp++ = *pivot;
    if (!pred(*pivot, lower)) // search left
      kd_range_query<J>(first, pivot, lower, upper, outp);
    if (pred(*pivot, upper)) // search right
      kd_range_query<J>(next(pivot), last, lower, upper, outp);
  } else {
    copy_if(first, last, outp, [&](const TupleType& x){
      return within(x, lower, upper);
    });
  }
  return;
}

template <size_t I,
          typename Iter,
          typename TupleType,
          typename OutIter>
void kd_rq_iters(Iter first, Iter last,
                 const TupleType& lower,
                 const TupleType& upper,
                 OutIter outp)
{
  if (distance(first, last) > 32) {
    auto pred = less_nth<I>();
    auto pivot = middle_of(first, last);
    constexpr auto J = next_dim<I, TupleType>;
    if (within(*pivot, lower, upper)) *outp++ = pivot;
    if (!pred(*pivot, lower)) // search left
      kd_rq_iters<J>(first, pivot, lower, upper, outp);
    if (pred(*pivot, upper)) // search right
      kd_rq_iters<J>(next(pivot), last, lower, upper, outp);
  } else {
    while (first != last) {
      if(within(*first, lower, upper)) *outp++ = first;
      ++first;
    }
  }
  return;
}

template <size_t I,
          typename Iter,
          typename TupleType,
          typename OutIter>
void kd_range_query(Iter first, Iter last,
                    const TupleType& center,
                    double radius,
                    OutIter outp)
{
  if (distance(first, last) > 32) {
    auto pred = less_radius_nth<I>();
    auto pivot = middle_of(first, last);
    constexpr auto J = next_dim<I, TupleType>;
    if (l2dist(*pivot, center) <= radius) *outp++ = *pivot;
    if (!pred(*pivot, center, -radius)) // search left
      kd_range_query<J>(first, pivot, center, radius, outp);
    if (!pred(center, *pivot, -radius)) // search right
      kd_range_query<J>(next(pivot), last, center, radius, outp);
  } else {
    copy_if(first, last, outp, [&](const TupleType& x){
      return l2dist(x, center) <= radius;
    });
  }
  return;
}

template <size_t I,
          typename Iter,
          typename TupleType,
          typename OutIter>
void kd_rq_iters(Iter first, Iter last,
                 const TupleType& center,
                 double radius,
                 OutIter outp)
{
  if (distance(first, last) > 32) {
    auto pred = less_radius_nth<I>();
    auto pivot = middle_of(first, last);
    constexpr auto J = next_dim<I, TupleType>;
    if (l2dist(*pivot, center) <= radius) *outp++ = pivot;
    if (!pred(*pivot, center, -radius)) // search left
      kd_rq_iters<J>(first, pivot, center, radius, outp);
    if (pred(*pivot, center, radius)) // search right
      kd_rq_iters<J>(next(pivot), last, center, radius, outp);
  } else {
    while (first != last) {
      if (l2dist(*first, center) <= radius) *outp++ = first;
      ++first;
    }
  }
  return;
}

template <typename Key, typename Iter>
struct less_key {
  bool operator()(const pair<Key, Iter>& lhs, const pair<Key, Iter>& rhs) {
    return lhs.first < rhs.first;
  }
};

template <typename Iter, typename Key = double>
struct n_best
{
  using qmem_t = pair<Key, Iter>;
  using qcomp_t = less_key<Key, Iter>;
  using qcont_t = vector<qmem_t>;
  size_t m_n;
  qcont_t m_q;
  n_best(size_t n): m_n(n), m_q() {
    m_q.reserve(n);
  }
  Key max_key() const
  {
    return m_q.size() < m_n ?
    numeric_limits<Key>::max() :
    m_q.front().first;
  }
  void add(Key dist, Iter iter)
  {
    if (m_q.size() < m_n) {
      m_q.emplace_back(dist, iter);
      if (m_q.size() == m_n)
        make_heap(begin(m_q), end(m_q), qcomp_t());
    } else if (dist < m_q.front().first) {
      pop_heap(begin(m_q), end(m_q), qcomp_t());
      m_q.back().first = dist; m_q.back().second = iter;
      push_heap(begin(m_q), end(m_q), qcomp_t());
    }
  }
  template <typename OutIter>
  void copy_to(OutIter outp)
  {
    if (m_q.size() < m_n) sort(begin(m_q), end(m_q), qcomp_t());
    else sort_heap(begin(m_q), end(m_q), qcomp_t());
    transform(begin(m_q), end(m_q), outp, [](const qmem_t& x){
      return *x.second;
    });
  }
  template <typename OutIter>
  void copy_iters_to(OutIter outp)
  {
    if (m_q.size() < m_n) sort(begin(m_q), end(m_q), qcomp_t());
    else sort_heap(begin(m_q), end(m_q), qcomp_t());
    transform(begin(m_q), end(m_q), outp, [](const qmem_t& x){
      return x.second;
    });
  }
  template <typename OutIter>
  void copy_dist_to(OutIter outp)
  {
    if (m_q.size() < m_n) sort(begin(m_q), end(m_q), qcomp_t());
    else sort_heap(begin(m_q), end(m_q), qcomp_t());
    copy(begin(m_q), end(m_q), outp);
  }
};

template <size_t I,
          typename Iter,
          typename TupleType,
          typename QType>
void knn(Iter first, Iter last,
         const TupleType& key,
         QType& Q)
{
  switch(distance(first, last)) {
  case 1 : Q.add(l2dist(*first, key), first);
  case 0 : return; } // switch end
  auto pivot = middle_of(first, last);
  Q.add(l2dist(*pivot, key), pivot);
  constexpr auto J = next_dim<I, TupleType>;
  if (equal_nth<I>()(*pivot, key)) {
    knn<J>(first, pivot, key, Q);
    knn<J>(next(pivot), last, key, Q);
  } else {
    auto search_left = less_nth<I>()(key, *pivot);
    if (search_left)
      knn<J>(first, pivot, key, Q);
    else
      knn<J>(next(pivot), last, key, Q);
    if (dist_nth<I>(key, *pivot) <= Q.max_key())
    {
      if (search_left)
        knn<J>(next(pivot), last, key, Q);
      else
        knn<J>(first, pivot, key, Q);
    }
  }
}

template <size_t I,
          typename Iter,
          typename TupleType,
          typename QType>
void knn(Iter first, Iter last,
         const TupleType& key,
         double p,
         QType& Q)
{
  switch(distance(first, last)) {
  case 1 : Q.add(pdist(*first, key, p), first);
  case 0 : return; } // switch end
  auto pivot = middle_of(first, last);
  Q.add(pdist(*pivot, key, p), pivot);
  constexpr auto J = next_dim<I, TupleType>;
  if (equal_nth<I>()(*pivot, key)) {
    knn<J>(first, pivot, key, Q);
    knn<J>(next(pivot), last, key, Q);
  } else {
    auto search_left = less_nth<I>()(key, *pivot);
    if (search_left)
      knn<J>(first, pivot, key, Q);
    else
      knn<J>(next(pivot), last, key, Q);
    if (dist_nth<I>(key, *pivot) <= Q.max_key())
    {
      if (search_left)
        knn<J>(next(pivot), last, key, Q);
      else
        knn<J>(first, pivot, key, Q);
    }
  }
}

} // namespace detail

namespace utils {

using detail::all_less;
using detail::none_less;
using detail::within;

using detail::middle_of;
using detail::median_part;

using detail::is_last;
using detail::next_dim;
using detail::iter_value_t;

using detail::kd_less;
using detail::less_nth;
using detail::equal_nth;
using detail::pred_nth;
using detail::dist_nth;
using detail::kd_compare;
using detail::make_kd_compare;

using detail::pdist;
using detail::l2dist;

using detail::n_best;

} // namespace utils

template <typename Iter>
void lex_sort(Iter first, Iter last)
{
  std::sort(first, last, utils::kd_less<0>());
}

template <typename Iter, typename Compare>
void lex_sort(Iter first, Iter last, const Compare& comp)
{
  std::sort(first, last, utils::make_kd_compare<0>(comp));
}

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

template <typename Iter, typename Compare>
void kd_sort_threaded(Iter first, Iter last, const Compare& comp)
{
  detail::kd_sort_threaded<0>(first, last, comp);
}

template <typename Iter>
bool kd_is_sorted(Iter first, Iter last)
{
  return detail::kd_is_sorted<0>(first, last);
}

template <typename Iter, typename Compare>
bool kd_is_sorted(Iter first, Iter last, const Compare& comp)
{
  return detail::kd_is_sorted<0>(first, last, comp);
}

template <typename Iter>
bool kd_is_sorted_threaded(Iter first, Iter last)
{
  return detail::kd_is_sorted_threaded<0>(first, last);
}

template <typename Iter, typename Compare>
bool kd_is_sorted_threaded(Iter first, Iter last, const Compare& comp)
{
  return detail::kd_is_sorted_threaded<0>(first, last, comp);
}

template <typename Iter, typename Key>
Iter kd_lower_bound(Iter first, Iter last, const Key& key)
{
  return detail::kd_lower_bound<0>(first, last, key);
}

template <typename Iter, typename Key>
Iter kd_upper_bound(Iter first, Iter last, const Key& key)
{
  return detail::kd_upper_bound<0>(first, last, key);
}

template <typename Iter, typename TupleType>
bool kd_binary_search(Iter first, Iter last, const TupleType& key)
{
  first = detail::kd_lower_bound<0>(first, last, key);
  return first != last && utils::none_less(key, *first);
}

template <typename Iter, typename Key>
std::pair<Iter, Iter> kd_equal_range(Iter first, Iter last, const Key& key)
{
  return std::make_pair(detail::kd_lower_bound<0>(first, last, key),
                        detail::kd_upper_bound<0>(first, last, key));
}

template <typename Iter, typename TupleType>
Iter kd_nearest_neighbor(Iter first, Iter last, const TupleType& key)
{
  return detail::kd_nearest_neighbor<0>(first, last, key);
}

template <typename Iter, typename TupleType>
Iter kd_nearest_neighbor(Iter first, Iter last,
                         const TupleType& key,
                         double p)
{
  return detail::kd_nearest_neighbor<0>(first, last, key, p);
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

template <typename Iter,
          typename TupleType,
          typename OutIter>
void kd_rq_iters(Iter first, Iter last,
                 const TupleType& lower,
                 const TupleType& upper,
                 OutIter outp)
{
  detail::kd_rq_iters<0>(first, last, lower, upper, outp);
}

template <typename Iter,
          typename TupleType,
          typename OutIter>
void kd_range_query(Iter first, Iter last,
                    const TupleType& center,
                    double radius, OutIter outp)
{
  detail::kd_range_query<0>(first, last, center, radius, outp);
}

template <typename Iter,
          typename TupleType,
          typename OutIter>
void kd_rq_iters(Iter first, Iter last,
                 const TupleType& center,
                 double radius, OutIter outp)
{
  detail::kd_rq_iters<0>(first, last, center, radius, outp);
}

template <typename Iter,
          typename TupleType,
          typename OutIter>
void kd_nearest_neighbors(Iter first, Iter last, const TupleType& key, size_t n, OutIter outp)
{
  size_t m = std::distance(first, last);
  detail::n_best<Iter> Q(std::min(n, m));
  detail::knn<0>(first, last, key, Q);
  Q.copy_to(outp);
}

template <typename Iter,
          typename TupleType,
          typename OutIter>
void kd_nn_iters(Iter first, Iter last, const TupleType& key, size_t n, OutIter outp)
{
  size_t m = std::distance(first, last);
  detail::n_best<Iter> Q(std::min(n, m));
  detail::knn<0>(first, last, key, Q);
  Q.copy_iters_to(outp);
}

template <typename Iter,
          typename TupleType,
          typename OutIter>
void kd_nn_dist(Iter first, Iter last, const TupleType& key, size_t n, OutIter outp)
{
  size_t m = std::distance(first, last);
  detail::n_best<Iter> Q(std::min(n, m));
  detail::knn<0>(first, last, key, Q);
  Q.copy_dist_to(outp);
}

template <typename Iter,
          typename TupleType,
          typename OutIter>
void kd_nearest_neighbors(Iter first, Iter last,
                          const TupleType& key,
                          double p, size_t n,
                          OutIter outp)
{
  size_t m = std::distance(first, last);
  detail::n_best<Iter> Q(std::min(n, m));
  detail::knn<0>(first, last, key, p, Q);
  Q.copy_to(outp);
}

template <typename Iter,
          typename TupleType,
          typename OutIter>
void kd_nn_iters(Iter first, Iter last,
                 const TupleType& key,
                 double p, size_t n,
                 OutIter outp)
{
  size_t m = std::distance(first, last);
  detail::n_best<Iter> Q(std::min(n, m));
  detail::knn<0>(first, last, key, p, Q);
  Q.copy_iters_to(outp);
}

template <typename Iter,
          typename TupleType,
          typename OutIter>
void kd_nn_dist(Iter first, Iter last,
                const TupleType& key,
                double p, size_t n,
                OutIter outp)
{
  size_t m = std::distance(first, last);
  detail::n_best<Iter> Q(std::min(n, m));
  detail::knn<0>(first, last, key, p, Q);
  Q.copy_dist_to(outp);
}

} // namespace kdtools
} // namespace keittlab

#endif // NO_CXX17

#endif // KDTOOLS_H
