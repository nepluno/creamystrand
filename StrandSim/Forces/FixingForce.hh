/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef FIXINGFORCE_HH_
#define FIXINGFORCE_HH_

#include "../Core/ElasticStrand.hh"
#include "ViscousOrNotViscous.hh"

namespace strandsim {

template <typename ViscousT = NonViscous>
class FixingForce : public ForceBase {
 public:
  FixingForce() {}
  virtual ~FixingForce() {}

 public:
  static const IndexType s_first =
      0;  // The first index on which this force can apply
  static const IndexType s_last = 0;  // The last index (counting from the end)

  typedef Eigen::Matrix<Scalar, 4, 1> LocalForceType;
  typedef Eigen::Matrix<Scalar, 4, 4> LocalJacobianType;
  typedef VecXx ForceVectorType;
  typedef Scalar LocalMultiplierType;

  static std::string getName() { return ViscousT::getName() + "fixing"; }

  static Scalar localEnergy(const ElasticStrand& strand, StrandState& geometry,
                            const IndexType vtx);

  static void computeLocal(LocalMultiplierType& localL,
                           const ElasticStrand& strand, const IndexType vtx);

  static void computeLocal(LocalForceType& localF, const ElasticStrand& strand,
                           StrandState& geometry, const IndexType vtx);

  static void computeLocal(LocalJacobianType& localJ,
                           const ElasticStrand& strand, StrandState& geometry,
                           const IndexType vtx);

  static void addInPosition(ForceVectorType& globalForce, const IndexType vtx,
                            const LocalForceType& localForce);

  static void addInPosition(JacobianMatrixType& globalJacobian,
                            const IndexType vtx,
                            const LocalJacobianType& localJacobian);

  static void addInPositionMultiplier(VecXx& globalMultiplier,
                                      const IndexType vtx,
                                      const LocalMultiplierType& localL);

  static void accumulateCurrentE(Scalar& energy, ElasticStrand& strand);
  static void accumulateCurrentF(VecXx& force, ElasticStrand& strand);
  static void accumulateCurrentJ(JacobianMatrixType& Jacobian,
                                 ElasticStrand& strand);
};

}  // namespace strandsim

#endif /* FixingForce_HH_ */
