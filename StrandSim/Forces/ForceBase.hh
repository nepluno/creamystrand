/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef FORCEBASE_HH_
#define FORCEBASE_HH_

#include <list>
#include <memory>

#include "../Core/Definitions.hh"

namespace strandsim {

class ElasticStrand;
class StrandState;

class ForceBase {
 public:
  enum Quantities { NONE, E = 1, F = 2, EF = 3, J = 4, EFJ = 7 };
};

//! Force that depend on the state of other strands
class CoupledForceBase : public ForceBase {
 public:
};

/**
 * This is the base class for run-time external forces that are stored in an
 * ElasticStrand's polymorphic list m_runtimeForces. Each accumulate* method
 * accumulates energy, force and/or Jacobian in geometry
 */
class RuntimeForceBase : public ForceBase {
 public:
  RuntimeForceBase();
  virtual ~RuntimeForceBase();

  virtual std::string getName() const = 0;

  virtual void accumulate(Quantities q, StrandState& geometry,
                          const ElasticStrand& strand);

  virtual void accumulateEFJ(StrandState& geometry,
                             const ElasticStrand& strand) = 0;

  virtual void accumulateEF(StrandState& geometry,
                            const ElasticStrand& strand) = 0;

  virtual void accumulateE(StrandState& geometry, const ElasticStrand& strand);

  virtual void accumulateJ(StrandState& geometry,
                           const ElasticStrand& strand) = 0;

  virtual void accumulateCurrentF(VecXx& globalF,
                                  const ElasticStrand& strand) = 0;

  virtual bool reverseHairdoFlag(bool isCombtool) const = 0;

  // New style interface (consistent with static-class forces)
  virtual void accumulateFutureE(Scalar& energy,
                                 const ElasticStrand& strand) = 0;

  virtual void draw() const {}  // Debug

  // Self managing forces will add and remove themselves from each strand force
  // list ( Instead of being added to all strands by an external StrandManager )
  // They will also have to take care of restarting affected strands when they
  // change Typical example is SpringForce, which only affect a few strands

 private:
  StrandSet* m_affectedStrands;
};

//! Force that adds/removes itself to/from strands
class SelfManagingRuntimeForceBase : public RuntimeForceBase {
 public:
  SelfManagingRuntimeForceBase();
  virtual ~SelfManagingRuntimeForceBase();

  virtual void addSelf(const std::vector<ElasticStrand*>& strands) = 0;

  virtual void removeSelf(const std::vector<ElasticStrand*>& strands) = 0;
};

}  // namespace strandsim

#endif /* FORCEBASE_HH_ */
