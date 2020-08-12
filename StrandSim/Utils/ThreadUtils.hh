/**
 * \copyright 2012 Gilles Daviet, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef THREADUTILS_HH_
#define THREADUTILS_HH_

#include <Eigen/Core>
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <numeric>
#include <thread>
#include <vector>

namespace strandsim {

typedef boost::mutex MutexType;
typedef boost::lock_guard<MutexType> LockGuard;
typedef boost::unique_lock<MutexType> UniqueLock;
typedef boost::condition_variable Condition;

typedef boost::thread ThreadType;

/*! Since boost::mutex are non copyable, this class allows other to have a
 member mutex without having to redefine their copy constructor and assignement
 operator
 */
class MutexWrapper {
 public:
  MutexWrapper() {}
  MutexWrapper(const MutexWrapper&) {}
  MutexWrapper& operator=(const MutexWrapper&) { return *this; }

  MutexType& operator*() { return m_mutex; }

 private:
  MutexType m_mutex;
};

class ThreadHandle {
 public:
  ThreadHandle() : m_running(false) {}

  ~ThreadHandle() { join(); }

  template <typename DataT, void (DataT::*func)(void)>
  void run(DataT* callee) {
    while (m_running) {
      join();
    }

    LockGuard lock(m_mutex);
    callable<DataT, func> f;

    m_thread = boost::thread(f, callee);
    m_running = true;
  }

  template <typename DataT>
  void run(DataT* callee) {
    run<DataT, &DataT::operator()>(callee);
  }

  void join() {
    LockGuard lock(m_mutex);

    if (m_running) {
      m_thread.join();
      m_running = false;
    }
  }

 private:
  bool m_running;
  boost::thread m_thread;
  MutexType m_mutex;

  template <typename DataT, void (DataT::*func)(void)>
  struct callable {
    void operator()(DataT* callee) { (callee->*func)(); }
  };
};

template <typename ContainerT, typename ReturnT, typename CallableT,
          typename ArgT>
void for_each(ContainerT& list, CallableT& obj,
              ReturnT (CallableT::*func)(ArgT)) {
  auto callee = std::bind1st(std::mem_fun(func), &obj);
  for_each(list.begin(), list.end(), callee);
}

template <typename ContainerT, typename ReturnT, typename CallableT,
          typename ArgT>
void for_each(ContainerT& list, const CallableT& obj,
              ReturnT (CallableT::*func)(ArgT) const) {
  auto callee = std::bind1st(std::mem_fun(func), &obj);
  for_each(list.begin(), list.end(), callee);
}

template <typename ContainerT, typename ReturnT, typename CallableT,
          typename ArgT>
void parfor(ContainerT& vec, CallableT& obj, ReturnT (CallableT::*func)(ArgT)) {
  const unsigned N = vec.size();
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL && defined(_OPENMP)
#pragma omp parallel for
#endif
  for (unsigned i = 0; i < N; ++i) {
    (obj.*func)(vec[i]);
  }
}

template <typename ContainerT, typename ReturnT, typename CallableT,
          typename ArgT>
void parfor(ContainerT& vec, const CallableT& obj,
            ReturnT (CallableT::*func)(ArgT) const) {
  const int N = (int)vec.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (unsigned i = 0; i < N; ++i) {
    (obj.*func)(vec[i]);
  }
}

template <typename ContainerT, typename CallableT>
void parfor(ContainerT& vec, CallableT& obj) {
  const int N = (int)vec.size();

#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL && defined(_OPENMP)
#pragma omp parallel for
#endif
  for (unsigned i = 0; i < N; ++i) {
    obj(vec[i]);
  }
}

inline unsigned get_num_threads() {
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL
  return std::thread::hardware_concurrency();
#else
  return 1U;
#endif
}

template <typename Index, typename Callable>
void for_each(Index start, Index end, Callable func) {
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL && defined(_OPENMP)
#pragma omp parallel for
#endif
  for (Index i = start; i < end; ++i) {
    func(i);
  }
}

template <typename Data, typename Callable>
void for_each(const std::vector<Data>& vec, Callable func) {
  const size_t N = vec.size();
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL && defined(_OPENMP)
#pragma omp parallel for
#endif
  for (size_t i = 0; i < N; ++i) {
    func(vec[i]);
  }
}

template <typename Data, typename Callable>
void for_each(std::vector<Data>& vec, Callable func) {
  const size_t N = vec.size();
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL && defined(_OPENMP)
#pragma omp parallel for
#endif
  for (size_t i = 0; i < N; ++i) {
    func(vec[i]);
  }
}

template <typename Index, typename Data, typename Callable>
Data reduction(Data* array, Index n, const Data& base, Callable func) {
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL && defined(_OPENMP)
  std::vector<Data> res_buf(omp_get_max_threads(), base);
#pragma omp parallel for
  for (Index i = 0; i < n; ++i) {
    res_buf[omp_get_thread_num()] =
        func(res_buf[omp_get_thread_num()], array[i]);
  }

  Data res = res_buf[0];
  const Index num_threads = (Index)omp_get_max_threads();
  for (Index i = 1; i < num_threads; ++i) {
    res = func(res, res_buf[i]);
  }
  return res;
#else
  Data res = base;
  for (Index i = 0; i < n; ++i) {
    res = func(res, array[i]);
  }
  return res;
#endif
}

template <typename T>
void parallel_concatenate(std::vector<T>& result,
                          const std::vector<std::vector<T> >& input) {
  if (!input.size()) return;

  std::vector<int> numvals(input.size());
  for_each(0, (int)input.size(),
           [&](int idx) { numvals[idx] = (int)input[idx].size(); });

  std::partial_sum(numvals.begin(), numvals.end(), numvals.begin());
  const int total_num_vals = numvals[numvals.size() - 1];

  result.resize(total_num_vals);
  if (total_num_vals == 0) return;

  for_each(0, (int)input.size(), [&](int idx) {
    const int base_idx = (idx == 0) ? 0 : numvals[idx - 1];
    const int local_num_vals = (int)input[idx].size();
    for (int i = 0; i < local_num_vals; ++i) {
      result[base_idx + i] = input[idx][i];
    }
  });
}

template <typename T, const int N>
void parallel_concatenate(
    Eigen::Matrix<T, Eigen::Dynamic, 1>& result,
    const std::vector<std::vector<Eigen::Matrix<T, N, 1> > >& input) {
  if (!input.size()) return;

  std::vector<int> numvals(input.size());
  for_each(0, (int)input.size(),
           [&](int idx) { numvals[idx] = (int)input[idx].size(); });

  std::partial_sum(numvals.begin(), numvals.end(), numvals.begin());
  const int total_num_vals = numvals[numvals.size() - 1];

  result.resize(total_num_vals * N);
  if (total_num_vals == 0) return;

  for_each(0, (int)input.size(), [&](int idx) {
    const int base_idx = (idx == 0) ? 0 : numvals[idx - 1];
    const int local_num_vals = (int)input[idx].size();
    for (int i = 0; i < local_num_vals; ++i) {
      result.template segment<N>((base_idx + i) * N) = input[idx][i];
    }
  });
}
}  // namespace strandsim

#endif /* THREADUTILS_HH_ */
