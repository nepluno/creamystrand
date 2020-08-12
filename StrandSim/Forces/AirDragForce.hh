/**
 * \copyright 2014 Danny Kaufman
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSIM_AIRDRAGFORCE_HH
#define STRANDSIM_AIRDRAGFORCE_HH

#include "../Core/BandMatrixFwd.hh"
#include "ForceBase.hh"

namespace strandsim {

class AirDragForce : public ForceBase {
 public:
  static const IndexType s_first =
      0;  // The first index on which this force can apply
  static const IndexType s_last = 0;  // The last index (counting from the end)

  typedef Vec3x LocalForceType;
  typedef Mat3x LocalJacobianType;
  typedef VecXx ForceVectorType;

  virtual std::string getName() const { return "air drag"; }

  static void computeLocal(LocalForceType& localF, const ElasticStrand& strand,
                           const StrandState& geometry, const IndexType vtx);

  static void computeLocal(LocalJacobianType& localJ,
                           const ElasticStrand& strand,
                           const StrandState& geometry, const IndexType vtx);

  static void addInPosition(ForceVectorType& globalForce, const IndexType vtx,
                            const LocalForceType& localForce);

  static void addInPosition(JacobianMatrixType& globalJacobian,
                            const IndexType vtx,
                            const LocalJacobianType& localJacobian);

  static void setFrameVelocities(const Vec3x& Omega, const Vec3x& velOrigin);

 protected:
  AirDragForce();
  virtual ~AirDragForce();

  static Vec3x s_velOrigin;
  static Mat3x s_Omega_Cross;
};

}  // namespace strandsim

#endif  // STRANDSIM_AIRDRAGFORCE_HH
