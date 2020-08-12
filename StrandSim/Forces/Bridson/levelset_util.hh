/**
 * \copyright 2008 Robert Bridson
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LEVELSET_UTIL_H
#define LEVELSET_UTIL_H

namespace strandsim {
namespace bridson {
double fraction_inside(double phi_left, double phi_right);
double fraction_inside(double phi_bl, double phi_br, double phi_tl,
                       double phi_tr);
}  // namespace bridson
}  // namespace strandsim

#endif
