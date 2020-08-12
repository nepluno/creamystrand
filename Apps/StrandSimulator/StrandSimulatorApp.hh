/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef __StrandSim__StandSimulatorApp__
#define __StrandSim__StandSimulatorApp__

#include <tclap/CmdLine.h>

#include <iostream>

#ifdef WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <sys/time.h>
#endif

#include "ProblemStepper.hh"

/** Converts from a string to the given type */
template <class T>
inline void fromString(T& t, const std::string& str) {
  std::stringstream ss(str);
  ss >> t;
}

/** Converts the given input to a string */
template <class T>
inline std::string toString(const T& t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

template <typename T>
class ProblemConstraint : public TCLAP::Constraint<T> {
 public:
  ProblemConstraint(const T& lower, const T& upper)
      : m_lower(lower), m_upper(upper) {
    m_typeDesc = toString(lower) + "-" + toString(upper);
  }

  virtual std::string description() const { return m_typeDesc; }

  virtual std::string shortID() const { return m_typeDesc; }

  virtual bool check(const T& value) const {
    return ((value >= m_lower) && (value <= m_upper));
  }

 protected:
  T m_lower;
  T m_upper;
  std::string m_typeDesc;
};

#endif /* __StrandSim__StandSimulatorApp__ */
