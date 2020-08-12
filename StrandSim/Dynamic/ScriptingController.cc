/**
 * \copyright 2010 Breannan Smith
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "ScriptingController.hh"

namespace strandsim {

ScriptingController::ScriptingController() : m_time(0.), m_dt(1.) {}

ScriptingController::ScriptingController(double time, double dt)
    : m_time(time), m_dt(dt) {}

ScriptingController::~ScriptingController() {}

void ScriptingController::setTime(double time) { m_time = time; }

void ScriptingController::setDt(double dt) { m_dt = dt; }

double ScriptingController::getTime() const { return m_time; }

double ScriptingController::getDt() const { return m_dt; }

}  // namespace strandsim
