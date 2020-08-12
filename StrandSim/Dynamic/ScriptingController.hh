/**
 * \copyright 2010 Breannan Smith
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSIM_SCRIPTINGCONTROLLER_HH
#define STRANDSIM_SCRIPTINGCONTROLLER_HH

#include "../Core/Definitions.hh"

namespace strandsim {

class ScriptingController {
 public:
  ScriptingController();
  ScriptingController(Scalar time, Scalar dt);

  virtual ~ScriptingController();

  //! Steps forward the ScriptingController to the current value of m_time
  virtual bool execute() = 0;

  void setTime(Scalar time);
  void setDt(Scalar dt);

  Scalar getTime() const;
  Scalar getDt() const;

 protected:
  Scalar m_time;
  Scalar m_dt;
};

}  // namespace strandsim

#endif  // STRANDSIM_SCRIPTINGCONTROLLER_HH
