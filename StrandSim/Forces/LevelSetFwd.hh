/**
 * \copyright 2012 Gilles Daviet
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LEVELSETFWD_HH_
#define LEVELSETFWD_HH_

namespace strandsim {

template <int deg>
class InterpolatedLevelSet;

// Change interpolation degree here
#define LEVELSET_INTERPOLATION_DEGREE 1
typedef InterpolatedLevelSet<LEVELSET_INTERPOLATION_DEGREE> LevelSet;

}  // namespace strandsim

#endif /* LEVELSETFWD_HH_ */
