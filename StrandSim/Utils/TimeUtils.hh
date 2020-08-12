/**
 * \copyright 2012 Gilles Daviet
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TIMEUTILS_HH_
#define TIMEUTILS_HH_

#ifndef WIN32
#include <sys/time.h>
#else
#define NOMINMAX
#include <Windows.h>
#endif
#include <time.h>

#include <iostream>

namespace strandsim {

template <typename StreamT = std::ostream>
class TimerBase {
 public:
  TimerBase(const std::string& name, StreamT* out = NULL,
            const char* id_str = NULL)
      : m_out(out), m_name(name), m_id(0), m_id_str(NULL) {
    restart(id_str);
  }

  ~TimerBase() { print(); }

  void setOutStream(StreamT* out) { m_out = out; }

#ifdef WIN32
  void convert_filetime(struct timeval* out_tv,
                        const FILETIME* filetime) const {
    // Microseconds between 1601-01-01 00:00:00 UTC and 1970-01-01 00:00:00 UTC
    static const uint64_t EPOCH_DIFFERENCE_MICROS = 11644473600000000ull;

    // First convert 100-ns intervals to microseconds, then adjust for the
    // epoch difference
    uint64_t total_us = (((uint64_t)filetime->dwHighDateTime << 32) |
                         (uint64_t)filetime->dwLowDateTime) /
                        10;
    total_us -= EPOCH_DIFFERENCE_MICROS;

    // Convert to (seconds, microseconds)
    out_tv->tv_sec = (time_t)(total_us / 1000000);
    out_tv->tv_usec = (long)(total_us % 1000000);
  }
#endif

  double elapsed() const {
    struct timeval cur;
#ifdef WIN32
    FILETIME ft;
    GetSystemTimeAsFileTime((LPFILETIME)&ft);
    convert_filetime(&cur, &ft);
#else
    gettimeofday(&cur, 0);
#endif

    return (cur.tv_sec - m_start.tv_sec) * 1.e3 +
           (cur.tv_usec - m_start.tv_usec) * 1.e-3;
  }

  void print() const {
    if (m_out) {
      print(*m_out);
      *m_out << "\n";
      m_out->flush();
    }
  }

  template <typename OutStream>
  void print(OutStream& out) const {
    const double el = elapsed();
    out << m_name << " [";

    if (m_id_str)
      out << m_id_str;
    else
      out << m_id;

    out << "] => " << el << " ms ";
  }

  void restart(const char* id_str = NULL) {
    if (m_id) print();
    ++m_id;
    m_id_str = id_str;
#ifdef WIN32
    FILETIME ft;
    GetSystemTimeAsFileTime((LPFILETIME)&ft);
    convert_filetime(&m_start, &ft);
#else
    gettimeofday(&m_start, 0);
#endif
  }

 private:
  StreamT* m_out;

  std::string m_name;
  unsigned m_id;
  const char* m_id_str;

  struct timeval m_start;
};

class Timer : public TimerBase<std::ostream> {
 public:
  Timer(const std::string& name, bool verbose = true, const char* id_str = NULL)
      : TimerBase<std::ostream>(name, verbose ? &std::cerr : NULL, id_str) {}
};

template <typename StreamT>
inline std::ostream& operator<<(std::ostream& out,
                                const TimerBase<StreamT>& timer) {
  timer.print(out);
  return out;
}

}  // namespace strandsim

#endif /* TIMEUTILS_HH_ */
