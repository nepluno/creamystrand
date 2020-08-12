/**
 * \copyright 2011 Mark Leone
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "MsgInfo.hh"

#include <iostream>

namespace strandsim {

/// Output severity and message id to a stream (with a trailing colon and
/// space).  If the severity and message id are empty, nothing is printed.
std::ostream& operator<<(std::ostream& out, const MsgInfo& info) {
  const char* severityStr = MsgInfo::SeverityToString(info.GetSeverity());
  if (*severityStr && !info.GetId().empty()) {
    out << severityStr << " (" << info.GetId() << "): ";
  } else if (*severityStr) {
    out << severityStr << ": ";
  } else if (!info.GetId().empty()) {
    out << "(" << info.GetId() << "): ";
  }
  return out;
}

}  // namespace strandsim
