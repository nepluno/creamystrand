/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "FluidPressureForce.hh"

#include "../Core/BandMatrix.hh"
#include "../Core/ElasticStrand.hh"
#include "../Core/ElasticStrandUtils.hh"
#include "../Dynamic/FluidScriptingController.hh"
#include "../Dynamic/ImplicitStepper.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Utils/MathUtilities.hh"

namespace strandsim {

std::shared_ptr<FluidScriptingController> FluidPressureForce::s_controller =
    NULL;

FluidPressureForce::FluidPressureForce() {}

FluidPressureForce::~FluidPressureForce() {}

void FluidPressureForce::computeLocal(LocalForceType& localF,
                                      const ElasticStrand& strand,
                                      StrandState& geometry,
                                      const IndexType vtx) {
  if (!s_controller || !strand.getStepper() ||
      !s_controller->isStrandSubmerged(strand.getGlobalIndex()) ||
      strand.isVertexFreezed(vtx)) {
    localF = LocalForceType::Zero();
    return;
  }

  const Vec3x x_s0 = strand.getVertex(vtx);

  const Scalar phi = s_controller->get_dense_liquid_signed_distance(x_s0);
  const Scalar dx = s_controller->getDx();
  if (phi > -0.5 * dx) {
    localF = LocalForceType::Zero();
    return;
  }

  Vec3x gp = s_controller->get_pressure_gradient(x_s0);

  const Scalar len = strand.getVoronoiLength(vtx);
  const Scalar radiusA = strand.getRadiusA(vtx);
  const Scalar radiusB = strand.getRadiusB(vtx);
  const Scalar A = strand.getFutureFlowDOFArea(vtx);
  const Scalar vol = (M_PI * radiusA * radiusB + A) * len;
  Vec3x force = -gp * vol;

  Scalar max_len = strand.collisionParameters().m_impulseMaxNorm;
  if (max_len > 0. && force.norm() > max_len) {
    force = force.normalized() * max_len;
  }

  localF = force;

  // check_isnan("FPF_localF", localF);
}

void FluidPressureForce::computeLocal(LocalJacobianType& localJ,
                                      const ElasticStrand& strand,
                                      StrandState& geometry,
                                      const IndexType vtx) {
  //        if(!s_controller || !strand.getStepper() ||
  //        !s_controller->isStrandSubmerged(strand.getGlobalIndex()) ||
  //        strand.isVertexFreezed(vtx)) {
  localJ = LocalJacobianType::Zero();
  //            return;
  //        }
  //        const Vec3x x_s0 = strand.getVertex( vtx );
  //////        const Scalar phi =
  ///s_controller->get_dense_liquid_signed_distance(x_s0);
  //////        const Scalar dx = s_controller->getDx();
  //////        if( phi > -2.0 * dx )
  //////        {
  //////            localJ = LocalJacobianType::Zero();
  //////            return;
  //////        }
  //////
  ////
  ////
  //        const Mat3x gp = s_controller->get_pressure_hessian(x_s0);
  //        const Scalar len = strand.getVoronoiLength( vtx );
  //        const Scalar radiusA = strand.getRadiusA( vtx );
  //        const Scalar radiusB = strand.getRadiusB( vtx );
  //        const Scalar vol = M_PI * radiusA * radiusB * len;
  //        localJ = -gp * vol;

  // check_isnan("FPF_localJ", localJ);
}

void FluidPressureForce::setScriptingController(
    const std::shared_ptr<FluidScriptingController>& controller) {
  s_controller = controller;
}

void FluidPressureForce::addInPosition(ForceVectorType& globalForce,
                                       const IndexType vtx,
                                       const LocalForceType& localForce) {
  globalForce.segment<3>(4 * vtx) += localForce;
}

void FluidPressureForce::addInPosition(JacobianMatrixType& globalJacobian,
                                       const IndexType vtx,
                                       const LocalJacobianType& localJacobian) {
  globalJacobian.localStencilAdd<3>(4 * vtx, localJacobian);
}

} /* namespace strandsim */
