/**
 * \copyright 2011 Mark Leone
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MSG_INFO_HH_S
#define MSG_INFO_HH_S

#include <cassert>
#include <cstdarg>
#include <iosfwd>
#include <string>

namespace strandsim {

class MsgInfo {
 public:
  /// Message ids are represented as strings.
  typedef std::string Id;

  /// Message severity.
  enum Severity {
    kSilent,
    kTrace,
    kDebug,
    kCopious,
    kContact,
    kInfo,
    kNotice,
    kWarning,
    kError,
    kSevere,
    kInternal,
    kForced,
  };

  /// Message frequency.
  enum Frequency {
    kOncePerId,
    kOncePerMsg,
    kAlways,
  };

  /// The default severity level.
  static const Severity kDefaultSeverity = kNotice;

  /// Convert severity level to string.
  inline static const char* SeverityToString(Severity severity) {
    switch (severity) {
      case kSilent:
        return "";
      case kTrace:
        return "Trace";
      case kDebug:
        return "Debug";
      case kCopious:
        return "Copious";
      case kContact:
        return "Contact";
      case kInfo:
        return "Info";
      case kNotice:
        return "Notice";
      case kWarning:
        return "Warning";
      case kError:
        return "Error";
      case kSevere:
        return "Severe";
      case kInternal:
        return "Internal error";
      case kForced:
      default:
        return "";
    }
  }

  /// Construct a message identifier.
  MsgInfo(Severity severity, const Id& id, Frequency frequency = kAlways)
      : m_id(id), m_severity(severity), m_frequency(frequency) {}

  /// Get the message id.
  const Id& GetId() const { return m_id; }

  /// Get message severity.
  Severity GetSeverity() const { return m_severity; }

  /// Set message severity.
  void SetSeverity(Severity severity) { m_severity = severity; }

  /// Get the message frequency.
  Frequency GetFrequency() const { return m_frequency; }

 private:
  Id m_id;                ///< The message identifier.
  Severity m_severity;    ///< Message severity.
  Frequency m_frequency;  ///< Message frequency
};

/// Output severity and message id to a stream (with a trailing colon and
/// space).  If the severity and message id are empty, nothing is printed.
std::ostream& operator<<(std::ostream& out, const MsgInfo& info);

template <MsgInfo::Severity severity>
class SeverityInfo : public MsgInfo {
 public:
  explicit SeverityInfo(const Id& id, Frequency frequency = kAlways)
      : MsgInfo(severity, id, frequency) {}
};

typedef SeverityInfo<MsgInfo::kSilent> SilentMsg;
typedef SeverityInfo<MsgInfo::kTrace> TraceMsg;
typedef SeverityInfo<MsgInfo::kDebug> DebugMsg;
typedef SeverityInfo<MsgInfo::kCopious> CopiousMsg;
typedef SeverityInfo<MsgInfo::kInfo> InfoMsg;
typedef SeverityInfo<MsgInfo::kNotice> NoticeMsg;
typedef SeverityInfo<MsgInfo::kWarning> WarningMsg;
typedef SeverityInfo<MsgInfo::kError> ErrorMsg;
typedef SeverityInfo<MsgInfo::kSevere> SevereMsg;
typedef SeverityInfo<MsgInfo::kInternal> InternalMsg;
typedef SeverityInfo<MsgInfo::kForced> ForcedMsg;

}  // namespace strandsim

#endif  // ndef MSG_INFO_HH
