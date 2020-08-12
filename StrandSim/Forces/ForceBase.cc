/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "ForceBase.hh"

#include "../Core/StrandState.hh"
#include "../Utils/TextLog.hh"

namespace strandsim {

RuntimeForceBase::RuntimeForceBase() : ForceBase() {}

RuntimeForceBase::~RuntimeForceBase() {}

SelfManagingRuntimeForceBase::SelfManagingRuntimeForceBase()
    : RuntimeForceBase() {}

SelfManagingRuntimeForceBase::~SelfManagingRuntimeForceBase() {}

void RuntimeForceBase::accumulate(Quantities q, StrandState& geometry,
                                  const ElasticStrand& strand) {
  switch (q) {
    case EFJ:
      accumulateEFJ(geometry, strand);
      break;
    case EF:
    case F:
      accumulateEF(geometry, strand);
      break;
    case E:
      accumulateE(geometry, strand);
      break;
    case J:
      accumulateJ(geometry, strand);
      break;
    case NONE:
      break;
  }
}

void RuntimeForceBase::accumulateE(StrandState& geometry,
                                   const ElasticStrand& strand) {
  WarningStream(g_log, "")
      << " accumulateE not implemented for " << getName()
      << " -> falling back to accumulateEF ( inefficient ) ";
  const VecXx currentF = geometry.m_totalForce;
  accumulateEF(geometry, strand);
  geometry.m_totalForce = currentF;
}

}  // namespace strandsim
