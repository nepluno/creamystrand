/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "FluidDragForce.hh"

#include "../Core/BandMatrix.hh"
#include "../Core/ElasticStrand.hh"
#include "../Core/ElasticStrandUtils.hh"
#include "../Dynamic/FluidScriptingController.hh"
#include "../Dynamic/ImplicitStepper.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Utils/MathUtilities.hh"

namespace strandsim {

std::shared_ptr<FluidScriptingController> FluidDragForce::s_controller = NULL;

FluidDragForce::FluidDragForce() {}

FluidDragForce::~FluidDragForce() {}

void FluidDragForce::computeLocal(LocalForceType& localF,
                                  const ElasticStrand& strand,
                                  StrandState& geometry, const IndexType vtx) {
  if (!s_controller || !strand.getStepper() || !s_controller->doUseDrag() ||
      !s_controller->isStrandSubmerged(strand.getGlobalIndex()) ||
      strand.isVertexFreezed(vtx)) {
    localF = LocalForceType::Zero();
    if (s_controller)
      s_controller->storeDragCoeff(strand.getGlobalIndex(), vtx, 0.0);
    return;
  }

  const Vec3x x_s0 = strand.getCurrentState().getVertex(vtx);
  const Scalar phi = s_controller->get_liquid_phi(x_s0);
  const Scalar dx = s_controller->getDx();
  if (phi > -s_controller->getDragInsulation() * dx) {
    localF = LocalForceType::Zero();
    s_controller->storeDragCoeff(strand.getGlobalIndex(), vtx, 0.0);
    return;
  }

  const Vec3x u_s0 = strand.getStepper()->velocities().segment<3>(vtx * 4);
  const Vec3x u_s =
      strand.dynamics().getDisplacement(vtx) / strand.getParameters().getDt();
  const Vec3x u_f0 = s_controller->get_velocity(strand.getGlobalIndex(), vtx);
  const Vec3x edge = strand.getCurrentState().getEdgeVector(vtx);

  const Scalar len = strand.m_VoronoiLengths[vtx];
  const Scalar radius =
      strand.getRadiusA(vtx);  // + strand.getCurrentFlowHeight( vtx );
  const VecXx color = s_controller->getComponentsAtPos(x_s0);

  const Scalar rho = s_controller->getDensity(color);
  const Scalar mu = s_controller->getViscosity(color);
  const Scalar n = s_controller->getFlowBehaviorIndex(color);
  const Scalar tau_y = s_controller->getYieldStress(color);
  const Scalar epsilon = s_controller->getEpsilon(x_s0);

  Scalar coeff =
      computeLocal(localF, edge, u_s0, u_s, u_f0, len, radius, rho, mu, n,
                   tau_y, epsilon, s_controller->useConstantDrag());
  // check_isnan("fdf_coeff", coeff);
  // check_isnan("fdf_local_F", localF);

  s_controller->storeDragCoeff(strand.getGlobalIndex(), vtx, coeff);
}

void FluidDragForce::computeLocal(LocalJacobianType& localJ,
                                  const ElasticStrand& strand,
                                  StrandState& geometry, const IndexType vtx) {
  if (!s_controller || !strand.getStepper() || !s_controller->doUseDrag() ||
      !s_controller->isStrandSubmerged(strand.getGlobalIndex()) ||
      strand.isVertexFreezed(vtx)) {
    localJ = LocalJacobianType::Zero();
    return;
  }

  const Vec3x x_s0 = strand.getCurrentState().getVertex(vtx);
  const Scalar phi = s_controller->get_liquid_phi(x_s0);
  const Scalar dx = s_controller->getDx();
  if (phi > -s_controller->getDragInsulation() * dx) {
    localJ = LocalJacobianType::Zero();
    return;
  }

  const Vec3x u_s0 = strand.getStepper()->velocities().segment<3>(vtx * 4);
  const Vec3x u_s =
      strand.dynamics().getDisplacement(vtx) / strand.getParameters().getDt();
  const VecXx color = s_controller->getComponentsAtPos(x_s0);

  const Vec3x u_f0 = s_controller->get_velocity(strand.getGlobalIndex(), vtx);
  const Vec3x edge = strand.getCurrentState().getEdgeVector(vtx);
  const Scalar len = strand.m_VoronoiLengths[vtx];
  const Scalar radius =
      strand.getRadiusA(vtx);  // + strand.getCurrentFlowHeight( vtx );
  const Scalar rho = s_controller->getDensity(color);
  const Scalar mu = s_controller->getViscosity(color);
  const Scalar n = s_controller->getFlowBehaviorIndex(color);
  const Scalar tau_y = s_controller->getYieldStress(color);
  const Scalar epsilon = s_controller->getEpsilon(x_s0);

  computeLocal(localJ, edge, u_s0, u_f0, len, radius, rho, mu, n, tau_y,
               epsilon, strand.getParameters().getDt(),
               s_controller->useConstantDrag());
}

void FluidDragForce::setScriptingController(
    const std::shared_ptr<FluidScriptingController>& controller) {
  s_controller = controller;
}

Scalar FluidDragForce::computeCoeff(
    const Vec3x& edge, const Vec3x& u_s0, const Vec3x& u_f0,
    const Scalar intersectLength, const Scalar radius, const Scalar rho,
    const Scalar viscosity, const Scalar behavior, const Scalar yield,
    const Scalar epsilon, bool useConstantDrag) {
  Vec3x dv = u_f0 - u_s0;
  const Scalar mdv = dv.norm();
  if (mdv == 0.0) {
    return 0.0;
  }

  const Scalar s2d3 = 0.8164965809;        // sqrt(2/3)
  const Scalar ln6 = 1.7917594692;         // ln(6)
  const Scalar ln3 = 1.0986122887;         // ln(3)
  const Scalar beta_coeff = 0.5613413994;  // 11/48*sqrt(6)
  const Scalar lns6m1ds6 = -0.5246681416;  // ln((sqrt(6) - 1)/sqrt(6))

  Vec3x ne = edge.normalized();
  const Scalar area_perp =
      intersectLength * radius * 2.0 * ne.cross(dv).norm() / mdv +
      M_PI * radius * radius;

  Scalar Cdi = 0.44;

  const Scalar dp = 2.0 * sqrt(area_perp / M_PI);  // |u-v|d^2, cm^3/s
  const Scalar Rei = std::max(
      rho * dp * epsilon * mdv * mdv /
          (viscosity * pow(mdv, behavior) + s2d3 * pow(dp, behavior) * yield),
      1e-2);
  const Scalar chi =
      3.7 - 0.65 * exp(-(1.5 - log10(Rei)) * (1.5 - log10(Rei)) * 0.5);

  if (!useConstantDrag) {
    const Scalar area_c =
        intersectLength * 2.0 * M_PI * radius + 4.0 * M_PI * radius * radius;
    const Scalar alpha = 3.0 / (behavior * behavior + behavior + 1.0);
    const Scalar X =
        pow(6.0, (behavior - 1.0) / 2.0) * pow(alpha, 1.0 + behavior);
    //        const Scalar Cdi = 0.44;
    //        const Scalar Cdi = 24.0 * X / Rei;
    Scalar Cd0 = 24.0 * X / Rei;

    const Scalar Cdinf = 0.44;
    const Scalar b = exp(3.0 * (alpha - ln6));
    const Scalar k = (3.0 - alpha) / (6.0 * alpha) *
                     exp((3.0 - alpha) / (2.0 * alpha) * ln3);
    const Scalar beta =
        beta_coeff * (1.0 - exp((3.0 - alpha) * (3.0 - alpha) /
                                (4.0 * alpha * alpha) * lns6m1ds6));

    Cdi = Cd0 +
          area_perp / area_c * Cdinf * pow(Cd0, 2.0 * beta) * k *
              pow((6.0 * X * b) / (6.0 * X * b + Cd0), beta) +
          Cdinf * (6.0 * X * b) / (6.0 * X * b + Cd0 * 128.0);
  }

  return Cdi * rho * area_perp * mdv * pow(epsilon, -chi) * 0.5;
}

Scalar FluidDragForce::computeLocal(
    LocalForceType& localF, const Vec3x& edge, const Vec3x& u_s0,
    const Vec3x& u_s, const Vec3x& u_f, const Scalar intersectLength,
    const Scalar radius, const Scalar rho, const Scalar viscosity,
    const Scalar behavior, const Scalar yield, const Scalar epsilon,
    bool useConstantDrag) {
  Vec3x dv = u_f - u_s;
  const Scalar mdv = dv.norm();
  Scalar coeff = 0.0;
  if (mdv == 0) {
    localF = Vec3x::Zero();
  } else {
    coeff = computeCoeff(edge, u_s0, u_f, intersectLength, radius, rho,
                         viscosity, behavior, yield, epsilon, useConstantDrag);
    localF = coeff * (u_f - u_s);
  }
  return coeff;
}

void FluidDragForce::computeLocal(LocalJacobianType& localJ, const Vec3x& edge,
                                  const Vec3x& u_s0, const Vec3x& u_f,
                                  const Scalar intersectLength,
                                  const Scalar radius, const Scalar rho,
                                  const Scalar viscosity, const Scalar behavior,
                                  const Scalar yield, const Scalar epsilon,
                                  const Scalar dt, bool useConstantDrag) {
  const Scalar coeff =
      computeCoeff(edge, u_s0, u_f, intersectLength, radius, rho, viscosity,
                   behavior, yield, epsilon, useConstantDrag);
  localJ = -coeff / dt * Mat3x::Identity();
}

void FluidDragForce::addInPosition(ForceVectorType& globalForce,
                                   const IndexType vtx,
                                   const LocalForceType& localForce) {
  globalForce.segment<3>(4 * vtx) += localForce;
}

void FluidDragForce::addInPosition(JacobianMatrixType& globalJacobian,
                                   const IndexType vtx,
                                   const LocalJacobianType& localJacobian) {
  globalJacobian.localStencilAdd<3>(4 * vtx, localJacobian);
}

} /* namespace strandsim */
