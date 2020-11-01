/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "BendingForce.hh"

#include "../Core/BandMatrix.hh"
#include "../Core/ElasticStrand.hh"
#include "../Core/ElasticStrandUtils.hh"
#include "../Dynamic/ImplicitStepper.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "ViscousOrNotViscous.hh"

namespace strandsim {

template <typename ViscousT>
void BendingForce<ViscousT>::computeLocal(
    typename BendingForce::LocalMultiplierType& localL,
    const ElasticStrand& strand, const IndexType vtx) {
  // B.ilen.(phi+hJv)
  const Mat2x& B = ViscousT::bendingMatrix(strand, vtx);
  const Vec4x& kappaBar = ViscousT::kappaBar(strand, vtx);
  const Scalar ilen = strand.m_invVoronoiLengths[vtx];
  const Vec4x& kappa = strand.m_currentState->m_kappas[vtx];
  const GradKType& gradKappa = strand.m_currentState->m_gradKappas[vtx];

  Vec4x Jvh = gradKappa.transpose() *
              strand.dynamics().getDisplacements().segment<11>(4 * (vtx - 1));

  localL.segment<2>(0) =
      -ilen * B *
      (kappa.segment<2>(0) - kappaBar.segment<2>(0) + Jvh.segment<2>(0));
  localL.segment<2>(2) =
      -ilen * B *
      (kappa.segment<2>(2) - kappaBar.segment<2>(2) + Jvh.segment<2>(2));
}

template <typename ViscousT>
Scalar BendingForce<ViscousT>::localEnergy(const ElasticStrand& strand,
                                           StrandState& geometry,
                                           const IndexType vtx) {
  const Mat2x& B = ViscousT::bendingMatrix(strand, vtx);
  const Vec4x& kappaBar = ViscousT::kappaBar(strand, vtx);
  const Scalar ilen = strand.m_invVoronoiLengths[vtx];
  const Vec4x& kappa = geometry.m_kappas[vtx];

  Scalar E =
      0.25 * ilen *
      ((kappa.segment<2>(0) - kappaBar.segment<2>(0))
           .dot(Vec2x(B * (kappa.segment<2>(0) - kappaBar.segment<2>(0)))) +
       (kappa.segment<2>(2) - kappaBar.segment<2>(2))
           .dot(Vec2x(B * (kappa.segment<2>(2) - kappaBar.segment<2>(2)))));

  // failsafes
  if (vtx > 0 && vtx < strand.getNumVertices() - 1) {
    const Scalar ks = ViscousT::kmb(strand, vtx) *
                      strand.getStepper()->getStretchMultiplier();
    const Scalar restLength = NonViscous::neighborDist(strand, vtx);

    const Scalar length = geometry.m_neighbor_distances[vtx];

    if (length < restLength) {
      const Scalar restState = ViscousT::neighborDist(strand, vtx);
      E += 0.5 * ks * square(length / restState - 1.) * restState;
    }
  }

  return E;
}

template <typename ViscousT>
void BendingForce<ViscousT>::computeLocal(
    typename BendingForce::LocalForceType& localF, const ElasticStrand& strand,
    StrandState& geometry, const IndexType vtx) {
  const Mat2x& B = ViscousT::bendingMatrix(strand, vtx);
  const Vec4x& kappaBar = ViscousT::kappaBar(strand, vtx);
  const Scalar ilen = strand.m_invVoronoiLengths[vtx];
  const Vec4x& kappa = geometry.m_kappas[vtx];
  const GradKType& gradKappa = geometry.m_gradKappas[vtx];

  localF = -ilen * 0.5 *
           (gradKappa.block<11, 2>(0, 0) * B *
                (kappa.segment<2>(0) - kappaBar.segment<2>(0)) +
            gradKappa.block<11, 2>(0, 2) * B *
                (kappa.segment<2>(2) - kappaBar.segment<2>(2)));

  if (vtx > 0 && vtx < strand.getNumVertices() - 1) {
    const Scalar l0 =
        (strand.getVertex(vtx + 1) - strand.getVertex(vtx - 1)).norm();
    const Scalar restLength = strand.m_restNeighborDistances[vtx];

    if (l0 < restLength) {
      const Scalar ks = ViscousT::kmb(strand, vtx) *
                        strand.getStepper()->getStretchMultiplier();
      const Scalar length = geometry.m_neighbor_distances[vtx];
      const Scalar restState = ViscousT::neighborDist(strand, vtx);

      const Vec3x edge = geometry.m_jump_edges[vtx] / length;
      Vec3x f = ks * (length / restState - 1.0) * edge;

      localF.segment<3>(0) += f;
      localF.segment<3>(8) += -f;
    }
  }
}

template <typename ViscousT>
void BendingForce<ViscousT>::computeLocal(
    typename BendingForce::LocalJacobianType& localJ,
    const ElasticStrand& strand, StrandState& geometry, const IndexType vtx) {
  localJ = geometry.m_bendingProducts[vtx] * 0.5;

  if (strand.m_requiresExactJacobian) {
    const Mat2x& bendingMatrixBase = strand.m_parameters.bendingMatrixBase();
    const Vec4x& kappaBar = ViscousT::kappaBar(strand, vtx);
    const Vec4x& kappa = geometry.m_kappas[vtx];
    const LocalJacobianType& hessKappa0 = geometry.m_hessKappas[vtx * 4 + 0];
    const LocalJacobianType& hessKappa1 = geometry.m_hessKappas[vtx * 4 + 1];
    const LocalJacobianType& hessKappa2 = geometry.m_hessKappas[vtx * 4 + 2];
    const LocalJacobianType& hessKappa3 = geometry.m_hessKappas[vtx * 4 + 3];
    const Vec2x& temp_e =
        bendingMatrixBase * (kappa.segment<2>(0) - kappaBar.segment<2>(0));
    const Vec2x& temp_f =
        bendingMatrixBase * (kappa.segment<2>(2) - kappaBar.segment<2>(2));

    localJ += (temp_e(0) * hessKappa0 + temp_e(1) * hessKappa1 +
               temp_f(0) * hessKappa2 + temp_f(1) * hessKappa3) *
              0.5;  
  }

  const Scalar ilen = strand.m_invVoronoiLengths[vtx];
  localJ *= -ilen * ViscousT::bendingCoefficient(strand, vtx);

  if (vtx > 0 && vtx < strand.getNumVertices() - 1) {
    const Scalar l0 =
        (strand.getVertex(vtx + 1) - strand.getVertex(vtx - 1)).norm();
    const Scalar restLength = strand.m_restNeighborDistances[vtx];

    if (l0 < restLength) {
      const Scalar ks = ViscousT::kmb(strand, vtx) *
                        strand.getStepper()->getStretchMultiplier();
      const Scalar length = geometry.m_neighbor_distances[vtx];
      const Scalar restState = ViscousT::neighborDist(strand, vtx);
      const Vec3x edge = geometry.m_jump_edges[vtx] / length;

      Mat3x M;

      if (strand.m_requiresExactJacobian) {
        M = ks * ((1.0 / restState - 1.0 / length) * Mat3x::Identity() +
                  1.0 / length * (edge * edge.transpose()));
      } else {
        M = ks *
            (std::max(0.0, 1.0 / restState - 1.0 / length) * Mat3x::Identity() +
             1.0 / length * (edge * edge.transpose()));
      }

      localJ.block<3, 3>(0, 0) += -M;
      localJ.block<3, 3>(8, 8) += -M;
      localJ.block<3, 3>(0, 8) += M;
      localJ.block<3, 3>(8, 0) += M;
    }
  }

  if (strand.m_projectJacobian) {
    Eigen::EigenSolver<typename BendingForce::LocalJacobianType> es; 
    es.compute(localJ, false);
    localJ -= std::max(0.0, es.eigenvalues().real().maxCoeff()) *
              LocalJacobianType::Identity();
  }
}

template <typename ViscousT>
void BendingForce<ViscousT>::addInPosition(VecXx& globalForce,
                                           const IndexType vtx,
                                           const LocalForceType& localForce) {
  globalForce.segment<11>(4 * (vtx - 1)) += localForce;
}

template <typename ViscousT>
void BendingForce<ViscousT>::addInPosition(
    JacobianMatrixType& globalJacobian, const IndexType vtx,
    const LocalJacobianType& localJacobian) {
  globalJacobian.localStencilAdd<11>(4 * (vtx - 1), localJacobian);
}

template <typename ViscousT>
void BendingForce<ViscousT>::accumulateCurrentE(Scalar& energy,
                                                ElasticStrand& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    energy += localEnergy(strand, strand.getCurrentState(), vtx);
  }
}

template <typename ViscousT>
void BendingForce<ViscousT>::accumulateCurrentF(VecXx& force,
                                                ElasticStrand& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    LocalForceType localF;
    computeLocal(localF, strand, strand.getCurrentState(), vtx);
    addInPosition(force, vtx, localF);
  }
}

template <typename ViscousT>
void BendingForce<ViscousT>::addInPositionMultiplier(
    VecXx& globalMultiplier, const IndexType vtx,
    const LocalMultiplierType& localL) {
  globalMultiplier.segment<4>(4 * vtx) += localL;
}

template <typename ViscousT>
void BendingForce<ViscousT>::accumulateCurrentJ(JacobianMatrixType& Jacobian,
                                                ElasticStrand& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    LocalJacobianType localJ;
    computeLocal(localJ, strand, strand.getCurrentState(), vtx);
    addInPosition(Jacobian, vtx, localJ);
  }
}

template class BendingForce<NonViscous>;
template class BendingForce<Viscous>;

}  // namespace strandsim
