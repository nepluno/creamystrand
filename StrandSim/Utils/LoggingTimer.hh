/**
 * \copyright 2011 Mark Leone
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LOGGINGTIMER_HH
#define LOGGINGTIMER_HH

#include "TextLog.hh"
#include "TimeUtils.hh"

namespace strandsim {

template <typename StreamName>
class LoggingTimer : public TimerBase<StreamCreator<StreamName> > {
 public:
  typedef StreamCreator<StreamName> StreamType;

  LoggingTimer(const std::string& name, const char* id_str = NULL)
      : TimerBase<StreamCreator<StreamName> >(name, &stream(), id_str) {}

  static StreamType& stream() {
    static StreamType s_stream(g_log, "Timing");
    return s_stream;
  }
};

}  // namespace strandsim

#endif  // LOGGINGTIMER_HH
