/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDBASE_H_
#define STRANDBASE_H_

#include "Definitions.hh"

namespace strandsim {

class StrandBaseParameters {};

class StrandBase {
  typedef StrandBaseParameters ParametersType;

 public:
  StrandBase();
  virtual ~StrandBase();
};

}  // namespace strandsim

#endif /* STRANDBASE_H_ */
