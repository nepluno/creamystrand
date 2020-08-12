/**
 * \copyright 2010 Breannan Smith
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "MeshScriptingController.hh"

namespace strandsim {

MeshScriptingController::MeshScriptingController(double time, double dt)
    : ScriptingController(time, dt) {}

MeshScriptingController::~MeshScriptingController() {}

bool MeshScriptingController::execute() { return execute(false); }

}  // namespace strandsim
