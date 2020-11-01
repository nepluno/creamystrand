/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "StretchingForce.hh"

#include <Eigen/Eigenvalues>

#include "../Core/BandMatrix.hh"
#include "../Core/BandMatrixFwd.hh"
#include "../Core/ElasticStrand.hh"
#include "../Dynamic/ImplicitStepper.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "ViscousOrNotViscous.hh"

namespace strandsim {

template <typename ViscousT>
void StretchingForce<ViscousT>::computeLocal(LocalMultiplierType& localL,
                                             const ElasticStrand& strand,
                                             const IndexType vtx) {
  const Scalar ks = ViscousT::ks(strand, vtx) *
                    strand.getStepper()->getStretchMultiplier();  // dyne
  const Scalar restLength = ViscousT::ellBar(strand, vtx);

  const Scalar length = strand.m_currentState->m_lengths[vtx];
  const Vec3x& edge = strand.m_currentState->m_tangents[vtx];

  localL =
      -ks *
      (length - ViscousT::ellBar(strand, vtx) +
       edge.dot(strand.dynamics().getDisplacements().segment<3>((vtx + 1) * 4) -
                strand.dynamics().getDisplacements().segment<3>(
                    vtx * 4)));  // dyne.cm
}

template <typename ViscousT>
Scalar StretchingForce<ViscousT>::localEnergy(const ElasticStrand& strand,
                                              StrandState& geometry,
                                              const IndexType vtx) {
  const Scalar ks =
      ViscousT::ks(strand, vtx) * strand.getStepper()->getStretchMultiplier();
  const Scalar restLength = ViscousT::ellBar(strand, vtx);

  const Scalar length = geometry.m_lengths[vtx];

  return 0.5 * ks * square(length / restLength - 1.0) * restLength;
}

template <typename ViscousT>
void StretchingForce<ViscousT>::addInPositionMultiplier(
    VecXx& globalMultiplier, const IndexType vtx,
    const LocalMultiplierType& localL) {
  globalMultiplier(vtx) += localL;
}

template <typename ViscousT>
void StretchingForce<ViscousT>::computeLocal(
    typename StretchingForce::LocalForceType& localF,
    const ElasticStrand& strand, StrandState& geometry, const IndexType vtx) {
  const Scalar ks =
      ViscousT::ks(strand, vtx) * strand.getStepper()->getStretchMultiplier();
  const Scalar restLength = ViscousT::ellBar(strand, vtx);

  const Scalar length = geometry.m_lengths[vtx];

  const Vec3x& edge = geometry.m_tangents[vtx];

  Vec3x f = ks * (length / restLength - 1.0) * edge;
  localF.segment<3>(0) = f;
  localF.segment<3>(3) = -f;
}

template <typename ViscousT>
void StretchingForce<ViscousT>::computeLocal(
    typename StretchingForce::LocalJacobianType& localJ,
    const ElasticStrand& strand, StrandState& geometry, const IndexType vtx) {
  const Scalar ks =
      ViscousT::ks(strand, vtx) * strand.getStepper()->getStretchMultiplier();
  const Scalar restLength = ViscousT::ellBar(strand, vtx);

  const Scalar length = geometry.m_lengths[vtx];
  const Vec3x& edge = geometry.m_tangents[vtx];

  Mat3x M;

  if (strand.m_requiresExactJacobian && !strand.m_projectJacobian) {
    M = ks * ((1.0 / restLength - 1.0 / length) * Mat3x::Identity() +
              1.0 / restLength * (edge * edge.transpose()));

  } else {
    M = ks *
        (std::max(0.0, 1.0 / restLength - 1.0 / length) * Mat3x::Identity() +
         1.0 / restLength * (edge * edge.transpose()));
  }

  localJ.block<3, 3>(0, 0) = localJ.block<3, 3>(3, 3) = -M;
  localJ.block<3, 3>(0, 3) = localJ.block<3, 3>(3, 0) = M;
}

template <typename ViscousT>
void StretchingForce<ViscousT>::addInPosition(
    typename StretchingForce::ForceVectorType& globalForce, const IndexType vtx,
    const LocalForceType& localForce) {
  globalForce.segment<3>(4 * vtx) += localForce.segment<3>(0);
  globalForce.segment<3>(4 * (vtx + 1)) += localForce.segment<3>(3);
}

template <typename ViscousT>
void StretchingForce<ViscousT>::addInPosition(
    JacobianMatrixType& globalJacobian, const IndexType vtx,
    const LocalJacobianType& localJacobian) {
  globalJacobian.edgeStencilAdd<6>(4 * vtx, localJacobian);
}

template <typename ViscousT>
void StretchingForce<ViscousT>::accumulateCurrentE(Scalar& energy,
                                                   ElasticStrand& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    energy += localEnergy(strand, strand.getCurrentState(), vtx);
  }
}

template <typename ViscousT>
void StretchingForce<ViscousT>::accumulateCurrentF(VecXx& force,
                                                   ElasticStrand& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    LocalForceType localF;
    computeLocal(localF, strand, strand.getCurrentState(), vtx);
    addInPosition(force, vtx, localF);
  }
}

template <typename ViscousT>
void StretchingForce<ViscousT>::accumulateCurrentJ(JacobianMatrixType& Jacobian,
                                                   ElasticStrand& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    LocalJacobianType localJ;
    computeLocal(localJ, strand, strand.getCurrentState(), vtx);
    addInPosition(Jacobian, vtx, localJ);
  }
}

template class StretchingForce<NonViscous>;
template class StretchingForce<Viscous>;

}  // namespace strandsim
