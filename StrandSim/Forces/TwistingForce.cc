/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "TwistingForce.hh"

#include "../Core/BandMatrix.hh"
#include "../Core/ElasticStrand.hh"
#include "../Dynamic/ImplicitStepper.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "ViscousOrNotViscous.hh"

namespace strandsim {

template <typename ViscousT>
void TwistingForce<ViscousT>::computeLocal(LocalMultiplierType& localL,
                                           const ElasticStrand& strand,
                                           const IndexType vtx) {
  const Scalar kt = ViscousT::kt(strand, vtx);
  const Scalar undefTwist = ViscousT::thetaBar(strand, vtx);
  const Scalar ilen = strand.m_invVoronoiLengths[vtx];
  const Scalar twist = strand.m_currentState->m_twists[vtx];

  localL =
      -kt * ilen *
      (twist - undefTwist +
       strand.m_currentState->m_gradTwists[vtx].dot(
           strand.dynamics().getDisplacements().segment<11>(4 * (vtx - 1))));
}

template <typename ViscousT>
Scalar TwistingForce<ViscousT>::localEnergy(const ElasticStrand& strand,
                                            StrandState& geometry,
                                            const IndexType vtx) {
  const Scalar kt =
      ViscousT::kt(strand, vtx) * strand.getStepper()->getStretchMultiplier();
  const Scalar undefTwist = ViscousT::thetaBar(strand, vtx);
  const Scalar ilen = strand.m_invVoronoiLengths[vtx];
  const Scalar twist = geometry.m_twists[vtx];

  return 0.5 * kt * square(twist - undefTwist) * ilen;
}

template <typename ViscousT>
void TwistingForce<ViscousT>::computeLocal(
    typename TwistingForce::LocalForceType& localF, const ElasticStrand& strand,
    StrandState& geometry, const IndexType vtx) {
  const Scalar kt =
      ViscousT::kt(strand, vtx) * strand.getStepper()->getStretchMultiplier();
  const Scalar undefTwist = ViscousT::thetaBar(strand, vtx);
  const Scalar ilen = strand.m_invVoronoiLengths[vtx];
  const Scalar twist = geometry.m_twists[vtx];

  localF = -kt * ilen * (twist - undefTwist) * geometry.m_gradTwists[vtx];
}

template <typename ViscousT>
void TwistingForce<ViscousT>::computeLocal(
    typename TwistingForce::LocalThetaForceType& localF,
    const ElasticStrand& strand, StrandState& geometry, const IndexType vtx) {
  const Scalar kt =
      ViscousT::kt(strand, vtx) * strand.getStepper()->getStretchMultiplier();
  const Scalar undefTwist = ViscousT::thetaBar(strand, vtx);
  const Scalar ilen = strand.m_invVoronoiLengths[vtx];
  const Scalar twist = geometry.m_twists[vtx];
  const Vec11x& gradTwist = geometry.m_gradTwists[vtx];
  Vec2x thetaGradTwist(gradTwist[3], gradTwist[7]);

  localF = -kt * ilen * (twist - undefTwist) * thetaGradTwist;
}

template <typename ViscousT>
void TwistingForce<ViscousT>::computeLocal(
    typename TwistingForce::LocalJacobianType& localJ,
    const ElasticStrand& strand, StrandState& geometry, const IndexType vtx) {
  const Scalar kt =
      ViscousT::kt(strand, vtx) * strand.getStepper()->getStretchMultiplier();
  const Scalar ilen = strand.m_invVoronoiLengths[vtx];
  const Mat11x& gradTwistSquared = geometry.m_gradTwistsSquared[vtx];

  localJ = -kt * ilen * gradTwistSquared;

  if (strand.m_requiresExactJacobian) {
    const Scalar undeformedTwist = ViscousT::thetaBar(strand, vtx);
    const Scalar twist = geometry.m_twists[vtx];
    const Mat11x& hessTwist = geometry.m_hessTwists[vtx];
    localJ += -kt * ilen * (twist - undeformedTwist) * hessTwist;
  }

  if (strand.m_projectJacobian) {
    Eigen::EigenSolver<typename TwistingForce::LocalJacobianType> es;
    es.compute(localJ, false);
    localJ -= std::max(0.0, es.eigenvalues().real().maxCoeff()) *
              LocalJacobianType::Identity();
  }
}

template <typename ViscousT>
void TwistingForce<ViscousT>::computeLocal(
    typename TwistingForce::LocalThetaJacobianType& localJ,
    const ElasticStrand& strand, StrandState& geometry, const IndexType vtx) {
  const Scalar kt =
      ViscousT::kt(strand, vtx) * strand.getStepper()->getStretchMultiplier();
  const Scalar undeformedTwist = ViscousT::thetaBar(strand, vtx);
  const Scalar ilen = strand.m_invVoronoiLengths[vtx];
  const Vec11x& gradTwist = geometry.m_gradTwists[vtx];
  Vec2x thetaGradTwist(gradTwist[3], gradTwist[7]);

  // There is no twist Hessian on the theta coordinate.
  localJ = -kt * ilen * (thetaGradTwist * thetaGradTwist.transpose());
}

template <typename ViscousT>
void TwistingForce<ViscousT>::addInPosition(VecXx& globalForce,
                                            const IndexType vtx,
                                            const LocalForceType& localForce) {
  globalForce.segment<11>(4 * (vtx - 1)) += localForce;
}

template <typename ViscousT>
void TwistingForce<ViscousT>::addInPosition(
    VecXx& globalForce, const IndexType vtx,
    const LocalThetaForceType& localForce) {
  globalForce.segment<2>(vtx - 1) += localForce;
}

template <typename ViscousT>
void TwistingForce<ViscousT>::addInPosition(
    JacobianMatrixType& globalJacobian, const IndexType vtx,
    const LocalJacobianType& localJacobian) {
  globalJacobian.localStencilAdd<11>(4 * (vtx - 1), localJacobian);
}

template <typename ViscousT>
void TwistingForce<ViscousT>::addInPosition(
    TriDiagonalMatrixType& globalJacobian, const IndexType vtx,
    const LocalThetaJacobianType& localJacobian) {
  globalJacobian.localStencilAdd<2>(vtx - 1, localJacobian);
}

template <typename ViscousT>
void TwistingForce<ViscousT>::accumulateCurrentE(Scalar& energy,
                                                 ElasticStrand& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    energy += localEnergy(strand, strand.getCurrentState(), vtx);
  }
}

template <typename ViscousT>
void TwistingForce<ViscousT>::accumulateCurrentF(VecXx& force,
                                                 ElasticStrand& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    LocalForceType localF;
    computeLocal(localF, strand, strand.getCurrentState(), vtx);
    addInPosition(force, vtx, localF);
  }
}

template <typename ViscousT>
void TwistingForce<ViscousT>::addInPositionMultiplier(
    VecXx& globalMultiplier, const IndexType vtx,
    const LocalMultiplierType& localL) {
  globalMultiplier(vtx) += localL;
}

template <typename ViscousT>
void TwistingForce<ViscousT>::accumulateCurrentJ(JacobianMatrixType& Jacobian,
                                                 ElasticStrand& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    LocalJacobianType localJ;
    computeLocal(localJ, strand, strand.getCurrentState(), vtx);
    addInPosition(Jacobian, vtx, localJ);
  }
}

template class TwistingForce<NonViscous>;
template class TwistingForce<Viscous>;

}  // namespace strandsim
