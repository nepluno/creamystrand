/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_TIMER_HPP
#define BOGUS_TIMER_HPP

#ifdef WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <sys/time.h>
#endif

namespace bogus {

//! Simple timer class. Starts when constructed.
class Timer {
 public:
  Timer() {
#ifdef WIN32
    LARGE_INTEGER proc_freq;
    ::QueryPerformanceFrequency(&proc_freq);
    m_freq = proc_freq.QuadPart;
#endif
    reset();
  }

  //! Resturns the elapsed time, in seconds, since the last call to reset()
  double elapsed() {
#ifdef WIN32
    LARGE_INTEGER stop;
    ::QueryPerformanceCounter(&stop);
    return ((stop.QuadPart - m_start.QuadPart) / m_freq);
#else
    struct timeval stop;
    gettimeofday(&stop, 0);
    return stop.tv_sec - m_start.tv_sec +
           1.e-6 * (stop.tv_usec - m_start.tv_usec);
#endif
  }

  //! Restarts the timer
  void reset() {
#ifdef WIN32
    ::QueryPerformanceCounter(&m_start);
#else
    gettimeofday(&m_start, 0);
#endif
  }

 private:
#ifdef WIN32
  double m_freq;
  LARGE_INTEGER m_start;
#else
  struct timeval m_start;
#endif
};

}  // namespace bogus

#endif
