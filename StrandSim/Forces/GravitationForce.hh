/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef GRAVITATIONFORCE_HH_
#define GRAVITATIONFORCE_HH_

#include "ForceBase.hh"

namespace strandsim {

class GravitationForce : public ForceBase {
 public:
  static const IndexType s_first =
      0;  // The first index on which this force can apply
  static const IndexType s_last = 0;  // The last index (counting from the end)

  typedef Vec3x LocalForceType;
  typedef Mat3x LocalJacobianType;

  GravitationForce();
  virtual ~GravitationForce();

  static std::string getName() { return "gravitation"; }

  static Scalar localEnergy(const ElasticStrand& strand,
                            const StrandState& geometry, const IndexType vtx);

  template <typename LocalT>
  static void computeLocal(LocalT& local, const ElasticStrand& strand,
                           const StrandState& geometry, const IndexType vtx);

  template <typename GlobalT, typename LocalT>
  static void addInPosition(GlobalT& global, const IndexType vtx,
                            const LocalT& local);

  static Vec3x getGravity() { return s_gravity; }

  static void setGravity(const Vec3x& gravity) { s_gravity = gravity; }

 private:
  static Vec3x s_gravity;
};

}  // namespace strandsim

#endif /* GRAVITATIONFORCE_HH_ */
