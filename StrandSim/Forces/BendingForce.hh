/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BENDINGFORCE_HH_
#define BENDINGFORCE_HH_

#include "ViscousOrNotViscous.hh"

namespace strandsim {

class ElasticStrand;

template <typename ViscousT = NonViscous>
class BendingForce : public ForceBase {
 public:
  BendingForce() {}
  virtual ~BendingForce() {}

 public:
  static const IndexType s_first =
      1;  // The first index on which this force can apply
  static const IndexType s_last = 1;  // The last index (counting from the end)

  typedef Eigen::Matrix<Scalar, 11, 1> LocalForceType;
  typedef Vec2x LocalThetaForceType;
  typedef Eigen::Matrix<Scalar, 11, 11> LocalJacobianType;
  typedef Mat2x LocalThetaJacobianType;
  typedef Eigen::Matrix<Scalar, 4, 1> LocalMultiplierType;

  static std::string getName() { return ViscousT::getName() + "bending"; }

  static Scalar localEnergy(const ElasticStrand& strand,
                            /*const*/ StrandState& geometry,
                            const IndexType vtx);

  static void computeLocal(LocalMultiplierType& localL,
                           const ElasticStrand& strand, const IndexType vtx);

  static void computeLocal(LocalForceType& localF, const ElasticStrand& strand,
                           StrandState& geometry, const IndexType vtx);
  static void computeLocal(LocalThetaForceType& localF,
                           const ElasticStrand& strand, StrandState& geometry,
                           const IndexType vtx);

  static void computeLocal(LocalJacobianType& localJ,
                           const ElasticStrand& strand, StrandState& geometry,
                           const IndexType vtx);

  static void addInPosition(VecXx& globalForce, const IndexType vtx,
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

#endif /* BENDINGFORCE_HH_ */
