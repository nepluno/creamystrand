/**
 * \copyright 2010 Breannan Smith, 2013 Danny Kaufman
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Timer.hh"

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

#ifdef WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <sys/time.h>
#endif

#include "../Core/Definitions.hh"

std::map<std::string, Timer*> Timer::s_timer_map;

// Static
Timer& Timer::getTimer(const std::string& name) {
  // If the timer already exists return it
  {
    std::map<std::string, Timer*>::iterator itr = s_timer_map.find(name);
    if (itr != s_timer_map.end()) {
      return *(itr->second);
    }
  }

  // Otherwise create a new timer and return it
  Timer* t = new Timer();
  s_timer_map[name] = t;
  return *t;
}

// Static
void Timer::deleteAllTimers() {
  for (std::map<std::string, Timer*>::iterator itr = s_timer_map.begin();
       itr != s_timer_map.end(); ++itr) {
    delete itr->second;
  }
  s_timer_map.clear();
}

// Static
void Timer::clearAllTimers() {
  for (std::map<std::string, Timer*>::iterator itr = s_timer_map.begin();
       itr != s_timer_map.end(); ++itr) {
    itr->second->clearTimer();
  }
}

// Static
void Timer::writeTimers(const std::string& file_name) {
  std::ofstream timer_stream(file_name.c_str(), std::ios::app);
  if (!timer_stream.is_open()) {
    std::cerr << "Error, failed to open timer log file: " << file_name
              << std::endl;
    return;
  }

  // Print the name of each timer
  for (std::map<std::string, Timer*>::iterator itr = s_timer_map.begin();
       itr != s_timer_map.end(); ++itr) {
    timer_stream << itr->first << "\t";
  }
  timer_stream << std::endl;

  // Print the value of each timer
  for (std::map<std::string, Timer*>::iterator itr = s_timer_map.begin();
       itr != s_timer_map.end(); ++itr) {
    timer_stream << itr->second->computeElapsedTime() << "\t";
  }
  timer_stream << std::endl;

  timer_stream.close();
}

Timer::Timer()
    : m_start_time(std::numeric_limits<strandsim::Scalar>::signaling_NaN()),
      m_elapsed_T(0),
      m_nesting_level(0),
      m_deepest_nesting_level(0) {}

void Timer::beginBlock() {
  assert(m_nesting_level >= 0);
  if (m_nesting_level == 0) {
    m_start_time = computeCurrentTime();
  }
  ++m_nesting_level;
  if (m_nesting_level > m_deepest_nesting_level) {
    m_deepest_nesting_level = m_nesting_level;
  }
}

void Timer::endBlock() {
  --m_nesting_level;
  assert(m_nesting_level >= 0);
  if (m_nesting_level == 0) {
    m_elapsed_T += (computeCurrentTime() - m_start_time);
  }
}

double Timer::computeElapsedTime() const {
  double elapsed_time = m_elapsed_T;

  // If active, add current elapsed time
  assert(m_nesting_level >= 0);
  if (m_nesting_level != 0) {
    elapsed_time += (computeCurrentTime() - m_start_time);
  }

  return elapsed_time;
}

double Timer::computeCurrentTime() const {
  timeval time;
#ifdef WIN32
  GetSystemTimeAsFileTime((LPFILETIME)&time);
#else
  gettimeofday(&time, NULL);
#endif
  return (double)time.tv_sec + ((double)time.tv_usec / 1000000.0);
}

void Timer::clearTimer() {
  m_start_time =
      std::numeric_limits<strandsim::Scalar>::signaling_NaN();  // SCALAR_NAN;
  m_elapsed_T = 0;
  m_nesting_level = 0;
  m_deepest_nesting_level = 0;
}
