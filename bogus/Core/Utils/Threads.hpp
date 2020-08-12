/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_UTILS_THREADS_HPP
#define BOGUS_UTILS_THREADS_HPP

#include "../Block/Constants.hpp"

#ifndef BOGUS_DONT_PARALLELIZE
#include <omp.h>
#endif

namespace bogus {

#ifdef BOGUS_DONT_PARALLELIZE
struct WithMaxThreads {
  explicit WithMaxThreads(int) {}
  int nThreads() const { return 1; }
};

#else
struct WithMaxThreads {
  explicit WithMaxThreads(int maxThreads)
      : m_prevMaxThreads(omp_get_max_threads()),
        m_newMaxThreads(maxThreads == 0 ? m_prevMaxThreads : maxThreads) {
    omp_set_num_threads(m_newMaxThreads);
  }

  ~WithMaxThreads() { omp_set_num_threads(m_prevMaxThreads); }

  int nThreads() const { return m_newMaxThreads; }

 private:
  const int m_prevMaxThreads;
  const int m_newMaxThreads;
};
#endif

}  // namespace bogus

#endif
