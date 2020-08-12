/**
 * \copyright 2011 Mark Leone
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TEXT_LOG_HH_S
#define TEXT_LOG_HH_S

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>

#include "MsgInfo.hh"
//#include <boost/thread/lock_guard.hpp>
#include <boost/timer.hpp>
#include <iosfwd>
#include <memory>
#include <unordered_map>
#include <unordered_set>

namespace strandsim {

/**
   TextLog provides message logging with configurable severity levels and
   duplicate suppression. Each message has a severity level and an identifer.
   Messages can be formatted via printf-style format strings:

       log->Write(ErrorMsg("E001"), "Expected value %f\n", x)

   Messages can also be formatted via streaming, allowing any value supporing
   the << operator to be written to a log:

       ErrorStream(log, "E001") << "Expected value " << x << "\n";

   By default, messages with low severity are suppressed.  The severity level of
   a message can be changed (e.g. to suppress a warning, or turn a warning into
   an error).  A message can also specify an optional frequency (see
   MsgInfo::Frequency), allowing duplicates to be suppressed either based on the
   message id or its full text.

   TextLog is thread-safe, and care is taken to avoid commingling the text of
   messages from different threads.
*/
class TextLog {
 public:
  /// Construct a log that writes to the given stream, suppressing messages
  /// below the specified level of severity.
  explicit TextLog(std::ostream& stream,
                   MsgInfo::Severity level = MsgInfo::kDefaultSeverity)
      : m_stream(stream), m_level(level), m_numErrors(0) {}

  /// Suppress messages below the specified level of severity.
  void SetMinSeverity(MsgInfo::Severity level) {
    boost::lock_guard<boost::mutex> lock(m_mutex);
    m_level = level;
  }

  /// Disable the message with the given id.
  void Disable(const MsgInfo::Id& id) { SetSeverity(id, MsgInfo::kSilent); }

  /// Override the severity of the message with the specified id.
  void SetSeverity(const MsgInfo::Id& id, MsgInfo::Severity severity) {
    boost::lock_guard<boost::mutex> lock(m_mutex);
    m_severities[id] = severity;
  }

  /// A log stream wraps a log in scoped lock for the duration of a single
  /// message and writes a newline at the end of the message.  Any value that
  /// supports operator<< can be written to a log stream.  A log stream can
  /// optionally discard all messages.
  /// TODO: Per-message suppresion for streams is not yet supported.
  class Stream {
   public:
    /// Construct a log stream from the given log for a message with the given
    /// info.  The stream might discard all output if the message is suppressed.
    /// Otherwise the log is locked via a scoped lock.
    Stream(TextLog* log, MsgInfo info) {
      // Lock the log until the destructor runs.
      m_mutex = log->getMutex();
      m_mutex->lock();

      // Note that the severity might be overridden.
      info.SetSeverity(log->getSeverity(info));
      m_out = log->getStream(info);
      if (m_out != NULL) *m_out << info;
    }

    /// The destructor writes a newline and unlocks the stream.
    ~Stream() {
      if (m_out) (*m_out) << std::endl;
      m_mutex->unlock();
    }

    /// Any value that supports operator<< can be written to a log.
    template <typename T>
    Stream& operator<<(const T& value) {
      if (m_out) {
        *m_out << value;
      }
      return *this;
    }

    /// Returns a non parallel-safe pointer to the stream
    static std::ostream* unsafe_stream(TextLog* log, MsgInfo info) {
      info.SetSeverity(log->getSeverity(info));
      std::ostream* out = log->getStream(info);
      if (out != NULL) *out << info;
      return out;
    }

    static bool would_print(TextLog* log, MsgInfo info) {
      return (log->getSeverity(info) >= log->m_level);
    }

   private:
    std::ostream* m_out;  // NULL if output is suppressed.
    boost::mutex* m_mutex;

    // Copying and assignment are prohibited.
    Stream(const Stream&);
    Stream& operator=(const Stream&);
  };

  /// Write a message using printf-style formatting, prefixed with its severity
  /// and message id, suffixed with a newline.  The message might be suppressed
  /// based on its severity and message id (see MsgInfo.hh).
  void Write(const MsgInfo& info, const char* format, ...) {
    va_list ap;
    va_start(ap, format);
    WriteV(info, format, ap);
    va_end(ap);
  }

  /// Write a message using varargs, prefixed with its severity and message id,
  /// suffixed with a newline.
  void WriteV(const MsgInfo& info, const char* format, va_list ap) {
    WriteWhereV(info, "", 0, format, ap);
  }

  /// Write a message that includes a filename and line number.  The filename is
  /// omitted if it is NULL or empty.  The line number is omitted if less or
  /// equal to zero.
  void WriteWhere(const MsgInfo& info, const char* filename, int line,
                  const char* format, ...) {
    va_list ap;
    va_start(ap, format);
    WriteWhereV(info, filename, line, format, ap);
    va_end(ap);
  }

  /// Varargs method for writing a message that includes a filename and line
  /// number.
  void WriteWhereV(MsgInfo info, const char* filename, int line,
                   const char* format, va_list ap);

  /// Write a string, without any locking or allocation (for signal handling
  /// safety).  The message is prefixed with the severity and message id and
  /// suffixed with a newline. Duplicate messages are not suppressed.
  void WriteSafely(const MsgInfo& info, const char* msg) const;

  /// Get the number of errors reported (excludes warnings, etc.)
  unsigned GetNumErrors() const { return m_numErrors; }

 protected:
  // Streams need access to protected methods.
  friend class Stream;

  /// Get the mutex.  All private methods require caller locking.
  boost::mutex* getMutex() const { return &m_mutex; }

  // Get the severity of the given message, which might have been overridden.
  MsgInfo::Severity getSeverity(const MsgInfo& info) {
    SeverityMap::const_iterator it = m_severities.find(info.GetId());
    return it == m_severities.end() ? info.GetSeverity() : it->second;
  }

  /// Get the output stream for a message with the given info.  Returns NULL if
  /// the message should be suppressed.
  std::ostream* getStream(const MsgInfo& info);

 private:
  // For now we format messages with vsnprintf using a fixed-sized local buffer.
  static const unsigned kMaxMsgSize = 4096;

  // Most operations lock the log for their duration.
  mutable boost::mutex m_mutex;

  // The output stream.
  std::ostream& m_stream;

  // Messages below this level of severity are suppressed.
  MsgInfo::Severity m_level;

  // Message severity can be overridden on a per-message basis (generalizing
  // -Woff, -Werror).
  typedef std::unordered_map<MsgInfo::Id, MsgInfo::Severity> SeverityMap;
  SeverityMap m_severities;

  // When requested, the message id is saved and subsequent messages with the
  // same id are suppressed.
  std::unordered_set<MsgInfo::Id> m_prevIds;

  // When requested, the message text is saved and subsequent identical messages
  // are suppressed.
  std::unordered_set<std::string> m_prevMsgs;

  // Number of errors reported (excludes warnings, etc.)
  unsigned m_numErrors;
};

/// This template provides a convenient syntax for constructing a log stream for
/// a message with a fixed severity.
template <MsgInfo::Severity severity>
class LogStream : public TextLog::Stream {
 public:
  LogStream(TextLog* log, const MsgInfo::Id& id,
            MsgInfo::Frequency frequency = MsgInfo::kAlways)
      : TextLog::Stream(log, MsgInfo(severity, id, frequency)) {}

  /// Returns a non parallel-safe pointer to the stream
  static std::ostream* unsafe_stream(
      TextLog* log, const MsgInfo::Id& id,
      MsgInfo::Frequency frequency = MsgInfo::kAlways) {
    return TextLog::Stream::unsafe_stream(log,
                                          MsgInfo(severity, id, frequency));
  }

  /// Returns whether a gioven message would be printed
  /// ( so we can avoid performing unnecessary computations )
  static bool would_print(TextLog* log, const MsgInfo::Id& id,
                          MsgInfo::Frequency frequency = MsgInfo::kAlways) {
    return TextLog::Stream::would_print(log, MsgInfo(severity, id, frequency));
  }
};

template <typename LogStreamT>
class StreamCreator {
 public:
  StreamCreator(TextLog* log, const MsgInfo::Id& id,
                MsgInfo::Frequency frequency = MsgInfo::kAlways)
      : m_log(log), m_id(id), m_frequency(frequency) {}

  ~StreamCreator() { flush(); }

  template <typename T>
  std::ostream& operator<<(const T& value) {
    static std::ostream blackhole(0);

    if (!m_stream) {
      // We have to use unsafe_stream because we do not wan't \n on destruction
      m_stream = LogStreamT::unsafe_stream(m_log, m_id, m_frequency);
    }
    if (m_stream) {
      return (*m_stream << value);
    }
    return blackhole;
  }

  void flush() { m_stream = NULL; }

 private:
  TextLog* m_log;
  MsgInfo::Id m_id;
  MsgInfo::Frequency m_frequency;

  std::ostream* m_stream;
};

typedef LogStream<MsgInfo::kSilent> SilentStream;
typedef LogStream<MsgInfo::kTrace> TraceStream;
typedef LogStream<MsgInfo::kDebug> DebugStream;
typedef LogStream<MsgInfo::kCopious> CopiousStream;
typedef LogStream<MsgInfo::kContact> ContactStream;
typedef LogStream<MsgInfo::kInfo> InfoStream;
typedef LogStream<MsgInfo::kNotice> NoticeStream;
typedef LogStream<MsgInfo::kWarning> WarningStream;
typedef LogStream<MsgInfo::kError> ErrorStream;
typedef LogStream<MsgInfo::kSevere> SevereStream;
typedef LogStream<MsgInfo::kInternal> InternalStream;
typedef LogStream<MsgInfo::kForced> ForcedStream;

extern TextLog* g_log;

}  // namespace strandsim

#endif  // ndef TEXT_LOG_HH
