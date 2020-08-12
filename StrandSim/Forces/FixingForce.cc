/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "FixingForce.hh"

#include "../Core/BandMatrix.hh"
#include "../Core/BandMatrixFwd.hh"
#include "../Core/ElasticStrand.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "ViscousOrNotViscous.hh"

namespace strandsim {

template <typename ViscousT>
void FixingForce<ViscousT>::computeLocal(LocalMultiplierType& localL,
                                         const ElasticStrand& strand,
                                         const IndexType vtx) {
  // not implemented!
}

template <typename ViscousT>
Scalar FixingForce<ViscousT>::localEnergy(const ElasticStrand& strand,
                                          StrandState& geometry,
                                          const IndexType vtx) {
  // not implemented!

  return 0.0;
}

template <typename ViscousT>
void FixingForce<ViscousT>::addInPositionMultiplier(
    VecXx& globalMultiplier, const IndexType vtx,
    const LocalMultiplierType& localL) {
  globalMultiplier(vtx) += localL;
}

template <typename ViscousT>
void FixingForce<ViscousT>::computeLocal(
    typename FixingForce::LocalForceType& localF, const ElasticStrand& strand,
    StrandState& geometry, const IndexType vtx) {
  localF.setZero();

  if (strand.isVertexGoaled(vtx)) {
    const Scalar ks = ViscousT::kf(strand, vtx);
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];

    if (ViscousT::isViscous()) {
      localF.segment<3>(0) =
          -ks * ilen * (geometry.getVertex(vtx) - strand.getVertex(vtx));
    } else {
      localF.segment<3>(0) =
          -ks * ilen * (geometry.getVertex(vtx) - strand.getGoalVertex(vtx));
    }
  }

  if (vtx < strand.getNumVertices() - 1 && strand.isThetaGoaled(vtx)) {
    const Scalar kt = ViscousT::kt(strand, vtx) *
                      strand.getParameters().getFixingMultiplier();
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];

    if (ViscousT::isViscous()) {
      localF(3) = -kt * ilen * (geometry.getTheta(vtx) - strand.getTheta(vtx));
    } else {
      localF(3) =
          -kt * ilen * (geometry.getTheta(vtx) - strand.getGoalTheta(vtx));
    }
  }
}

template <typename ViscousT>
void FixingForce<ViscousT>::computeLocal(
    typename FixingForce::LocalJacobianType& localJ,
    const ElasticStrand& strand, StrandState& geometry, const IndexType vtx) {
  localJ.setZero();

  if (strand.isVertexGoaled(vtx)) {
    const Scalar ks = ViscousT::kf(strand, vtx);
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];

    localJ.block<3, 3>(0, 0) = -ks * ilen * Mat3x::Identity();
  }

  if (vtx < strand.getNumVertices() - 1 && strand.isThetaGoaled(vtx)) {
    const Scalar kt = ViscousT::kt(strand, vtx) *
                      strand.getParameters().getFixingMultiplier();
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];

    localJ(3, 3) = -kt * ilen;
  }
}

template <typename ViscousT>
void FixingForce<ViscousT>::addInPosition(
    typename FixingForce::ForceVectorType& globalForce, const IndexType vtx,
    const LocalForceType& localForce) {
  const int num_verts = (globalForce.size() + 3) / 4;

  if (vtx < num_verts - 1) {
    globalForce.segment<4>(4 * vtx) += localForce;
  } else {
    globalForce.segment<3>(4 * vtx) += localForce.segment<3>(0);
  }
}

template <typename ViscousT>
void FixingForce<ViscousT>::addInPosition(
    JacobianMatrixType& globalJacobian, const IndexType vtx,
    const LocalJacobianType& localJacobian) {
  const int num_verts = (globalJacobian.rows() + 3) / 4;

  if (vtx < num_verts - 1) {
    globalJacobian.localStencilAdd<4>(4 * vtx, localJacobian);
  } else {
    Mat3x J = localJacobian.block<3, 3>(0, 0);
    globalJacobian.localStencilAdd<3>(4 * vtx, J);
  }
}

template <typename ViscousT>
void FixingForce<ViscousT>::accumulateCurrentE(Scalar& energy,
                                               ElasticStrand& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    energy += localEnergy(strand, strand.getCurrentState(), vtx);
  }
}

template <typename ViscousT>
void FixingForce<ViscousT>::accumulateCurrentF(VecXx& force,
                                               ElasticStrand& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    LocalForceType localF;
    computeLocal(localF, strand, strand.getCurrentState(), vtx);
    addInPosition(force, vtx, localF);
  }
}

template <typename ViscousT>
void FixingForce<ViscousT>::accumulateCurrentJ(JacobianMatrixType& Jacobian,
                                               ElasticStrand& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    LocalJacobianType localJ;
    computeLocal(localJ, strand, strand.getCurrentState(), vtx);
    addInPosition(Jacobian, vtx, localJ);
  }
}

template class FixingForce<NonViscous>;
template class FixingForce<Viscous>;

}  // namespace strandsim
