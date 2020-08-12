/**
 * \copyright 2011 Mark Leone
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "TextLog.hh"

#include <boost/thread/locks.hpp>
#include <cstdarg>
#include <cstdio>
#include <iostream>

namespace strandsim {

// Get the output stream for a message with the given info.  Returns NULL if the
// message should be suppressed.
std::ostream* TextLog::getStream(const MsgInfo& info) {
  // Suppress if the severity doesn't meet the specified level.
  MsgInfo::Severity severity = info.GetSeverity();
  if (severity < m_level) return NULL;

  // Supressed multiple occurrences of the same message id if requested.
  if (info.GetFrequency() == MsgInfo::kOncePerId) {
    bool isNew = m_prevIds.insert(info.GetId()).second;
    if (!isNew) return NULL;
  }

  // Keep a count of errors reported.
  if (severity == MsgInfo::kError || severity == MsgInfo::kSevere ||
      severity == MsgInfo::kInternal)
    ++m_numErrors;

  return &m_stream;
}

// Write a message, optionally suppressing duplicates.  The message is prefixed
// with the severity.
/// Varargs method for writing a message that includes a filename and line
/// number.
void TextLog::WriteWhereV(MsgInfo info, const char* filename, int line,
                          const char* format, va_list ap) {
  // Check whether the message should be suppressed.  If not, a lock is acquired
  // (the message info is not written).
  boost::mutex::scoped_lock lock;

  // Update the message severity, which might have been overridden.
  info.SetSeverity(getSeverity(info));

  // Get the output stream, or NULL if the message should be suppressed.
  std::ostream* stream = getStream(info);
  if (stream == NULL) return;

  // If we need to check for duplicate messages, write to a stringstream.
  std::ostream* out = stream;
  std::stringstream stringStream;
  bool suppressDuplicate = info.GetFrequency() == MsgInfo::kOncePerMsg;
  if (suppressDuplicate) out = &stringStream;

  // Write the filename and line number (if any)
  bool gotFilename = filename && *filename != '\0';
  if (gotFilename && line > 0)
    *out << filename << ':' << line << ": ";
  else if (gotFilename)
    *out << filename << ": ";
  else if (line > 0)
    *out << "line " << line << ": ";

  // Write the severity and message id (if any)
  *out << info;

  // Format and write the message.  For now we use vsnprintf with a fixed-sized
  // local buffer.
  // TODO: release lock while formatting?
  char msg[kMaxMsgSize];
  vsnprintf(msg, kMaxMsgSize, format, ap);
  *out << msg << std::endl;

  // Suppress identical duplicate messages if requested.
  if (suppressDuplicate) {
    bool isNew = m_prevMsgs.insert(stringStream.str()).second;
    if (!isNew) return;
    *stream << stringStream.str();
  }
}

// Write a string, without any locking or allocation (for signal handling
// safety).  The message is prefixed with the severity.
void TextLog::WriteSafely(const MsgInfo& info, const char* msg) const {
  // TODO: actually, it's not safe to write to a std::ostream from a signal
  // handler (nor flush). We could use a file-descriptor-based stream
  // implementation: http://www.josuttis.com/cppcode/fdstream.html
  if (info.GetSeverity() >= m_level) {
    m_stream << info << msg << std::endl;
  }
}

// Global variable. Must be initialized somewhere before any logging.
TextLog* g_log = NULL;

}  // namespace strandsim
