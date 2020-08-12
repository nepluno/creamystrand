/**
 * \copyright 2010 Breannan Smith, 2013 Danny Kaufman
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef __TIMER_H__
#define __TIMER_H__

#include <map>
#include <string>

// Macros to enable/disable timing at compile time

#ifdef TIMERS_ENABLED

#define START_TIMER(name) Timer::getTimer(name).beginBlock()

#define STOP_TIMER(name) Timer::getTimer(name).endBlock()

#define FETCH_TIMER(name) Timer::getTimer(name).computeElapsedTime()

#define SAVE_TIMERS(filename) Timer::writeTimers(filename)

#define DELETE_TIMERS() Timer::deleteAllTimers()

#define CLEAR_TIMERS() Timer::clearAllTimers()

#else

#define START_TIMER(name)

#define STOP_TIMER(name)

#define SAVE_TIMERS(filename)

#define DELETE_TIMERS()

#define CLEAR_TIMERS()

#endif

class Timer {
 public:
  Timer();

  // Returns the timer for a given name. Reuses the timer if it exists,
  // otherwise creates a new one.
  static Timer& getTimer(const std::string& name);

  // Deletes all timers created via the getTimer() function
  static void deleteAllTimers();

  // Sets all timers to 0
  static void clearAllTimers();

  // Outputs all timers to a file
  static void writeTimers(const std::string& file_name);

  // Start the timer by increasing the nesting level
  void beginBlock();

  // Decreases the nesting level, and stops the timer if the outermost block is
  // reached. When timer is stopped, the elapsed time is added the total elapsed
  // time.
  void endBlock();

  // Returns the total time
  double computeElapsedTime() const;

  void clearTimer();

 private:
  double computeCurrentTime() const;

  static std::map<std::string, Timer*> s_timer_map;

  // The last time the timer's start function was called
  double m_start_time;
  // Total time this timer has been active if not currently active.
  // Amount of time in previous matched begin/end calls if active.
  double m_elapsed_T;

  // If m_nesting_level is greater than zero, the timer is on
  int m_nesting_level;
  // The deepest nesting level encountered (helps spot unbalanced
  // begin/end-block)
  int m_deepest_nesting_level;
};

#endif
