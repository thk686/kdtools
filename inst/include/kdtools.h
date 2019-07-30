// Copyright Timothy H. Keitt 2018
// Please see the LICENSE file in the R package
// By using this header, you agree to cite this
// work in your product using the DOI found at
// https://doi.org/10.5281/zenodo.1213797

#ifndef __KDTOOLS_H__
#define __KDTOOLS_H__

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

namespace kdtools {

template <typename T>
struct ndim
{
  using U = typename std::remove_pointer<T>::type;
  static constexpr auto value = std::tuple_size<U>::value;
};

// Specialize for non-numeric types
// TODO: User defined distance functions

template <typename T>
double scalar_diff(const T& lhs, const T& rhs)
{
  return lhs - rhs;
}

template <typename T>
double scalar_dist(const T& lhs, const T& rhs)
{
  return std::abs(scalar_diff(lhs, rhs));
}

namespace detail {

using std::abs;
using std::get;
using std::next;
using std::prev;
using std::pair;
using std::size_t;
using std::thread;
using std::vector;
using std::is_same;
using std::distance;
using std::enable_if;
using std::partition;
using std::is_pointer;
using std::nth_element;
using std::tuple_element;
using std::numeric_limits;
using std::priority_queue;
using std::is_partitioned;
using std::remove_pointer;
using std::iterator_traits;
using std::partition_point;
using std::placeholders::_1;

template <typename T>
struct is_not_pointer
{
  static constexpr bool value = !is_pointer<T>::value;
};

template <typename T>
T middle_of(const T first, const T last)
{
  return next(first, distance(first, last) / 2);
}

template <size_t I>
struct less_nth
{
  template <typename T>
  typename enable_if<is_not_pointer<T>::value, bool>::type
  operator()(const T& lhs, const T& rhs)
  {
    return get<I>(lhs) < get<I>(rhs);
  }
  template <typename T>
  typename enable_if<is_pointer<T>::value, bool>::type
  operator()(const T& lhs, const T& rhs)
  {
    return get<I>(*lhs) < get<I>(*rhs);
  }
};

template <size_t I>
struct equal_nth
{
  template <typename T>
  typename enable_if<is_not_pointer<T>::value, bool>::type
  operator()(const T& lhs, const T& rhs)
  {
    return get<I>(lhs) == get<I>(rhs);
  }
  template <typename T>
  typename enable_if<is_pointer<T>::value, bool>::type
  operator()(const T& lhs, const T& rhs)
  {
    return get<I>(*lhs) == get<I>(*rhs);
  }
};

template <typename Pred, size_t I>
struct pred_nth
{
  Pred m_pred;
  pred_nth(const Pred& pred) : m_pred(pred) {}
  template <typename T>
  typename enable_if<is_not_pointer<T>::value, bool>::type
  operator()(const T& lhs, const T& rhs)
  {
    return m_pred(get<I>(lhs), get<I>(rhs));
  }
  template <typename T>
  typename enable_if<is_pointer<T>::value, bool>::type
  operator()(const T& lhs, const T& rhs)
  {
    return m_pred(get<I>(*lhs), get<I>(*rhs));
  }
};

template <size_t I, typename Pred>
pred_nth<Pred, I> make_pred_nth(const Pred& pred)
{
  return pred_nth<Pred, I>(pred);
}

template <size_t I, typename T>
typename enable_if<is_not_pointer<T>::value, double>::type
dist_nth(const T& lhs, const T& rhs)
{
  return scalar_dist(get<I>(lhs), get<I>(rhs));
}

template <size_t I, typename T>
typename enable_if<is_pointer<T>::value, double>::type
dist_nth(const T& lhs, const T& rhs)
{
  return scalar_dist(get<I>(*lhs), get<I>(*rhs));
}

template <size_t I, typename T>
typename enable_if<is_not_pointer<T>::value, double>::type
diff_nth(const T& lhs, const T& rhs)
{
  return scalar_diff(get<I>(lhs), get<I>(rhs));
}

template <size_t I, typename T>
typename enable_if<is_pointer<T>::value, double>::type
diff_nth(const T& lhs, const T& rhs)
{
  return scalar_diff(get<I>(*lhs), get<I>(*rhs));
}

template <size_t I, size_t N>
struct incr_wrap
{
  static constexpr auto value = (I + 1) % N;
};

template <size_t I, typename T>
struct next_dim
{
  static constexpr auto value = incr_wrap<I, ndim<T>::value>::value;
};

template <size_t I, typename TupleType>
struct is_not_last
{
  static constexpr auto value = (I != (ndim<TupleType>::value - 1));
};

template<size_t I, typename TupleType>
struct is_last
{
  static constexpr auto value = (I == (ndim<TupleType>::value - 1));
};

template <size_t I, size_t K = 0>
struct kd_less
{
  template <typename T>
  typename enable_if<is_not_last<K, T>::value, bool>::type
  operator()(const T& lhs, const T& rhs) const
  {
    constexpr auto J = next_dim<I, T>::value;
    return equal_nth<I>()(lhs, rhs) ?
      kd_less<J, K + 1>()(lhs, rhs) :
        less_nth<I>()(lhs, rhs);
  }
  template <typename T>
  typename enable_if<is_last<K, T>::value, bool>::type
  operator()(const T& lhs, const T& rhs) const
  {
    return less_nth<I>()(lhs, rhs);
  }
};

template <typename Pred, size_t I, size_t K = 0>
struct kd_compare
{
  Pred m_pred;
  kd_compare(const Pred& pred) : m_pred(pred) {}
  template <typename T>
  typename enable_if<is_not_last<K, T>::value, bool>::type
  operator()(const T& lhs, const T& rhs) const
  {
    auto pred = make_pred_nth<I>(m_pred);
    constexpr auto J = next_dim<I, T>::value;
    return !pred(lhs, rhs) && !pred(rhs, lhs) ?
      kd_compare<Pred, J, K + 1>(m_pred)(lhs, rhs) : pred(lhs, rhs);
  }
  template <typename T>
  typename enable_if<is_last<K, T>::value, bool>::type
  operator()(const T& lhs, const T& rhs) const
  {
    return make_pred_nth<I>(m_pred)(lhs, rhs);
  }
};

template <size_t I, typename Pred>
kd_compare<Pred, I> make_kd_compare(const Pred& pred)
{
  return kd_compare<Pred, I>(pred);
}

template <typename T>
using iter_value_t = typename iterator_traits<T>::value_type;

template <size_t I, typename Iter>
Iter find_pivot(Iter first, Iter last)
{
  using T = iter_value_t<Iter>;
  auto pivot = middle_of(first, last);
  return partition_point(first, pivot, [&](const T& x){
    return less_nth<I>()(x, *pivot);
  });
}

template <typename Iter, typename Pred>
Iter adjust_pivot(Iter first, Iter pivot, Pred pred)
{
  using T = iter_value_t<Iter>;
  return partition(first, pivot, [&](const T& x){
    return pred(x, *pivot);
  });
}

template <typename Iter, typename Pred>
Iter median_part(Iter first, Iter last, Pred pred)
{
  auto pivot = middle_of(first, last);
  nth_element(first, pivot, last, pred);
  return adjust_pivot(first, pivot, pred);
}

template <size_t I, typename Iter>
void kd_sort(Iter first, Iter last)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>::value;
  if (distance(first, last) > 1)
  {
    auto pred = kd_less<I>();
    auto pivot = median_part(first, last, pred);
    kd_sort<J>(next(pivot), last);
    kd_sort<J>(first, pivot);
  }
}

template <typename Iter, typename Pred>
bool check_partition(Iter first, Iter pivot, Iter last, Pred pred)
{
  using T = iter_value_t<Iter>;
  return is_partitioned(first, last, [&](const T& x){
    return pred(x, *pivot);
  });
}

template <size_t I, typename Iter>
bool kd_is_sorted(Iter first, Iter last)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>::value;
  if (distance(first, last) < 2) return true;
  auto pred = kd_less<I>();
  auto pivot = find_pivot<I>(first, last);
  return check_partition(first, pivot, last, pred) &&
    kd_is_sorted<J>(first, pivot) &&
    kd_is_sorted<J>(next(pivot), last);
}

template <size_t I, typename Iter, typename Compare>
void kd_sort(Iter first, Iter last, const Compare& comp)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>::value;
  if (distance(first, last) > 1)
  {
    auto pred = make_kd_compare<I>(comp);
    auto pivot = median_part(first, last, pred);
    kd_sort<J>(next(pivot), last, comp);
    kd_sort<J>(first, pivot, comp);
  }
}

template <size_t I, typename Iter, typename Compare>
bool kd_is_sorted(Iter first, Iter last, const Compare& comp)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>::value;
  if (distance(first, last) < 2) return true;
  auto pred = make_kd_compare<I>(comp);
  auto pivot = find_pivot<I>(first, last);
  return check_partition(first, pivot, last, pred) &&
    kd_is_sorted<J>(first, pivot, comp) &&
    kd_is_sorted<J>(next(pivot), last, comp);
}

template <size_t I, typename Iter>
void kd_sort_threaded(Iter first, Iter last,
                      int max_threads = std::thread::hardware_concurrency(),
                      int thread_depth = 1)
{
  using TupleType = iter_value_t<Iter>;
  constexpr auto J = next_dim<I, TupleType>::value;
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
  constexpr auto J = next_dim<I, TupleType>::value;
  if (distance(first, last) > 1)
  {
    auto pred = make_kd_compare<I>(comp);
    auto pivot = median_part(first, last, pred);
    if ((1 << thread_depth) <= max_threads)
    {
      thread t(kd_sort_threaded<J, Iter>,
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

template <size_t I>
struct all_less_
{
  template <typename T>
  typename enable_if<is_not_last<I, T>::value, bool>::type
  operator()(const T& lhs, const T& rhs) const
  {
    return less_nth<I>()(lhs, rhs) && all_less_<I + 1>()(lhs, rhs);
  }
  template <typename T>
  typename enable_if<is_last<I, T>::value, bool>::type
  operator()(const T& lhs, const T& rhs) const
  {
    return less_nth<I>()(lhs, rhs);
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
  typename enable_if<is_not_last<I, T>::value, bool>::type
  operator()(const T& lhs, const T& rhs)
  {
    return !less_nth<I>()(lhs, rhs) && none_less_<I + 1>()(lhs, rhs);
  }
  template <typename T>
  typename enable_if<is_last<I, T>::value, bool>::type
  operator()(const T& lhs, const T& rhs)
  {
    return !less_nth<I>()(lhs, rhs);
  }
};

template <typename TupleType>
bool none_less(const TupleType& lhs, const TupleType& rhs)
{
  return none_less_<0>()(lhs, rhs);
}

template <size_t I, typename Iter, typename TupleType>
Iter kd_lower_bound(Iter first, Iter last, const TupleType& value)
{
  constexpr auto J = next_dim<I, TupleType>::value;
  if (distance(first, last) > 1)
  {
    auto pivot = find_pivot<I>(first, last);
    if (none_less(*pivot, value))
      return kd_lower_bound<J>(first, pivot, value);
    if (all_less(*pivot, value))
      return kd_lower_bound<J>(next(pivot), last, value);
    auto it = kd_lower_bound<J>(first, pivot, value);
    if (it != last && none_less(*it, value)) return it;
    it = kd_lower_bound<J>(next(pivot), last, value);
    if (it != last && none_less(*it, value)) return it;
    return last;
  }
  return first != last && none_less(*first, value) ? first : last;
}

template <size_t I, typename Iter, typename TupleType>
Iter kd_upper_bound(Iter first, Iter last, const TupleType& value)
{
  constexpr auto J = next_dim<I, TupleType>::value;
  if (distance(first, last) > 1)
  {
    auto pivot = find_pivot<I>(first, last);
    if (all_less(value, *pivot))
      return kd_upper_bound<J>(first, pivot, value);
    if (none_less(value, *pivot))
      return kd_upper_bound<J>(next(pivot), last, value);
    auto it = kd_upper_bound<J>(first, pivot, value);
    if (it != last && all_less(value, *it)) return it;
    it = kd_upper_bound<J>(next(pivot), last, value);
    if (it != last && all_less(value, *it)) return it;
    return last;
  }
  return first != last && all_less(value, *first) ? first : last;
}

template <size_t I>
struct sum_of_squares_
{
  template <typename TupleType>
  typename enable_if<is_not_last<I, TupleType>::value, double>::type
  operator()(const TupleType& lhs, const TupleType& rhs) const
  {
    using next_ = sum_of_squares_<I + 1>;
    return std::pow(diff_nth<I>(rhs, lhs), 2) + next_()(lhs, rhs);
  }
  template <typename TupleType>
  typename enable_if<is_last<I, TupleType>::value, double>::type
  operator()(const TupleType& lhs, const TupleType& rhs) const
  {
    return std::pow(diff_nth<I>(rhs, lhs), 2);
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
  constexpr auto J = next_dim<I, TupleType>::value;
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
    if (dist_nth<I>(value, *pivot) < min_dist)
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
bool within(const TupleType& value,
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
  if (distance(first, last) > 32) {
    auto pred = less_nth<I>();
    auto pivot = find_pivot<I>(first, last);
    constexpr auto J = next_dim<I, TupleType>::value;
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
    auto pivot = find_pivot<I>(first, last);
    constexpr auto J = next_dim<I, TupleType>::value;
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

template <typename Iter, typename Key = double>
struct n_best
{
  using qmem_t = pair<Key, Iter>;
  using qcont_t = vector<qmem_t>;
  using qcomp_t = less_nth<0>;
  using queue_t = priority_queue<qmem_t, qcont_t, qcomp_t>;
  size_t m_n;
  queue_t m_q;
  n_best(size_t n) : m_n(n), m_q(qcomp_t()) {}
  Key max_key() const
  {
    return m_q.empty() ?
      numeric_limits<Key>::max() :
        m_q.top().first;
  }
  void add(Key dist, Iter it)
  {
    m_q.emplace(dist, it);
    if (m_q.size() > m_n) m_q.pop();
  }
  template <typename OutIter>
  void copy_to(OutIter outp)
  {
    while (!m_q.empty())
    {
      *outp++ = *m_q.top().second;
      m_q.pop();
    }
  }
  template <typename OutIter>
  void copy_iters_to(OutIter outp)
  {
    while (!m_q.empty())
    {
      *outp++ = m_q.top().second;
      m_q.pop();
    }
  }
};

template <size_t I,
          typename Iter,
          typename TupleType,
          typename QType>
void knn(Iter first, Iter last,
         const TupleType& value,
         QType& Q)
{
  switch(distance(first, last)) {
  case 1 : Q.add(l2dist(*first, value), first);
  case 0 : return; } // switch end
  auto pivot = find_pivot<I>(first, last);
  Q.add(l2dist(*pivot, value), pivot);
  auto search_left = less_nth<I>()(value, *pivot);
  constexpr auto J = next_dim<I, TupleType>::value;
  if (search_left)
    knn<J>(first, pivot, value, Q);
  else
    knn<J>(next(pivot), last, value, Q);
  if (dist_nth<I>(value, *pivot) <= Q.max_key())
  {
    if (search_left)
      knn<J>(next(pivot), last, value, Q);
    else
      knn<J>(first, pivot, value, Q);
  }
}

} // namespace detail

namespace utils {

using detail::all_less;
using detail::none_less;
using detail::within;

using detail::middle_of;
using detail::find_pivot;
using detail::median_part;

using detail::is_last;
using detail::is_not_last;
using detail::incr_wrap;
using detail::next_dim;

using detail::kd_less;
using detail::less_nth;
using detail::equal_nth;
using detail::pred_nth;
using detail::dist_nth;
using detail::kd_compare;
using detail::make_kd_compare;

using detail::l2dist;
using detail::sum_of_squares;

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

template <typename Iter>
bool kd_is_sorted(Iter first, Iter last)
{
  return detail::kd_is_sorted<0>(first, last);
}

template <typename Iter, typename Compare>
void kd_sort(Iter first, Iter last, const Compare& comp)
{
  detail::kd_sort<0>(first, last, comp);
}

template <typename Iter, typename Compare>
bool kd_is_sorted(Iter first, Iter last, const Compare& comp)
{
  return detail::kd_is_sorted<0>(first, last, comp);
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
  return first != last && utils::none_less(value, *first);
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
void kd_nearest_neighbors(Iter first, Iter last,
                          const TupleType& value,
                          size_t n, OutIter outp)
{
  detail::n_best<Iter> Q(n);
  detail::knn<0>(first, last, value, Q);
  Q.copy_to(outp);
}

template <typename Iter,
          typename TupleType,
          typename OutIter>
void kd_nn_iters(Iter first, Iter last,
                 const TupleType& value,
                 size_t n, OutIter outp)
{
  detail::n_best<Iter> Q(n);
  detail::knn<0>(first, last, value, Q);
  Q.copy_iters_to(outp);
}

} // namespace kdtools

#endif // __KDTOOLS_H__
