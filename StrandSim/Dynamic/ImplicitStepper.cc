/**
 * \copyright 2014 Danny Kaufman, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "ImplicitStepper.hh"

#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#include "../Core/ElasticStrand.hh"
#include "../Forces/ForceAccumulator.hh"
#include "../Forces/GravitationForce.hh"
#include "../Forces/StretchingForce.hh"
#include "../Utils/MathUtilities.hh"
#include "../Utils/TextLog.hh"
#include "DOFScriptingController.hh"
#include "FluidScriptingController.hh"
#include "SimulationParameters.hh"
#include "StrandDynamicTraits.hh"
#include "StrandImplicitManager.hh"
/*
 A dynamic step can trigger a succession of failsafes.
 He re is a summary of the process

 - startSubstep()

 Before contacts / hard constraints : solveUnconstrained() and update( false )
 - First try with linearized dynamics and exact jacobian.
 - If the jacobian is not SPD or the resulting step makes the strand stretch,
 use the non-exact jacobian, which sould be better conditionned.
 - If the jacobian is still not SPD or the strand is still stretching,
 fall back to the unconstraint non-linear solver ( solveNonLinear() )
 - If the strand is still badly stretching, use a geometric projection

 Optionally for contacts/hard-constraints:
 - Call prepareForExternalSolve()
 - Then either update( true ) to accept the step or rewind() to discard it
 - If update( true ) returns false, which mean the strand is stretching, the
 StramdImplicitManager will try other failsafe to re-solve the contacts and
 constraints

 - finalize()

 */

namespace strandsim {

ImplicitStepper::ImplicitStepper(ElasticStrand& strand,
                                 const SimulationParameters& params)
    : m_dt(0.),
      m_fraction(1.0),
      m_params(params),
      m_notSPD(false),
      m_usedNonlinearSolver(false),
      m_velocities(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
      m_savedVelocities(
          VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
      m_newVelocities(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
      m_beginningVelocities(
          VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
      m_flow_velocities(VecXx::Zero(strand.getNumEdges())),
      m_flow_newVelocities(VecXx::Zero(strand.getNumEdges())),
      m_flow_strain(VecXx::Zero(strand.getNumVertices())),
      m_flow_newStrain(VecXx::Zero(strand.getNumVertices())),
      m_projectionDisplacements(VecXx::Zero(m_velocities.rows())),
      m_additional_impulse(VecXx::Zero(m_velocities.rows())),
      m_additional_inertia(VecXx::Zero(m_velocities.rows())),
      m_strand(strand),
      m_newtonIter(0),
      m_total_flow(0.0),
      m_stretchMultiplier(1.0),
      m_stretchingFailureThreshold(10.0),
      m_costretchResidualFailureThreshold(1e-4) {
  const Scalar slenderness =
      m_strand.getTotalRestLength() /
      (m_strand.getRadiusA(0) *
       m_strand.getRadiusB(m_strand.getNumVertices() - 1));

  if (m_params.m_useImpulseMethod) {
    // set up mass matrix for zeroth order solve
    const unsigned ndofs = m_strand.getCurrentDegreesOfFreedom().rows();
    m_massMatrix.resize(ndofs, ndofs);
    m_massMatrix.setZero();
  }

  int num_components = 1;
  if (strand.getParent() && strand.getParent()->getFluidScriptingController()) {
    num_components =
        strand.getParent()->getFluidScriptingController()->getNumComponents();
  }

  const int num_verts = strand.getNumVertices();
  m_flow_components = VecXx::Zero(num_verts * num_components);

  // use zero component by default
  for (int i = 0; i < num_verts; ++i) {
    m_flow_components(i * num_components) = 1.0;
  }

  m_flow_newComponents = m_flow_components;
}

ImplicitStepper::~ImplicitStepper() {
  // for ( unsigned i = 0; i < m_bilateralConstraints.size(); ++i )
  // {
  //     delete m_bilateralConstraints[i];
  // }
  // delete m_nonlinearCallback ;
}

Scalar ImplicitStepper::getCurrentFlowVelocityAtVertex(int vtx) const {
  const int n = m_strand.getNumVertices();
  if (vtx == 0) {
    return m_flow_velocities[0];
  } else if (vtx == n - 1) {
    return m_flow_velocities[vtx - 1];
  } else {
    return (m_flow_velocities[vtx - 1] + m_flow_velocities[vtx]) * 0.5;
  }
}

void ImplicitStepper::updateDofController(int subStepId, int numSubstepsPerDt) {
  m_strand.dynamics().getScriptingController()->update(
      m_strand.getCurrentState(), subStepId, numSubstepsPerDt);
}

Scalar ImplicitStepper::getStretchMultiplier() const {
  return m_stretchMultiplier;
}

void ImplicitStepper::initStepping(Scalar dt) {
  this->m_dt = dt;
  // init everything at the beginning of a step that doesn't need substepping
  m_strand.setFutureAreaDegreesOfFreedom(
      m_strand.getCurrentAreaDegreesOfFreedom());

  // Save current state for substepping
  m_strand.setSavedDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom());

  // Only required if we read from checkpoint at previous t, but no easy way to
  // check that
  m_velocities = m_strand.dynamics().getDisplacements() / m_dt;
  m_beginningVelocities = m_velocities;
  m_savedVelocities = m_velocities;

  // Clear collisions from previous step
  m_additional_impulse.setZero();

  //		// check_isnan("vels", m_velocities);
}

void ImplicitStepper::startSubstep(int id, Scalar dt, Scalar fraction) {
  // Here we init everything may change during substepping
  this->m_dt = dt;
  this->m_fraction = fraction;

  m_notSPD = false;

  m_strand.requireExactJacobian(true);

  m_usedNonlinearSolver = false;

  m_strand.dynamics().computeViscousForceCoefficients(
      m_dt);  // Maybe an overkill to do this each time but at least we can
              // change the stepper's time step without worrying about it.
  m_strand.dynamics().getDisplacements() = m_velocities * m_dt;

  // recover stretch multiplier so that we can use default value as soon as
  // possible
  m_stretchMultiplier = std::min(1.0, m_stretchMultiplier * 1.5);
}

void ImplicitStepper::recoverFromBeginning() {
  // this must be executed after update()
  rewind();

  m_strand.setCurrentDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom());
  m_velocities = m_beginningVelocities;

  // refresh displacement and swap back.
  m_strand.dynamics().acceptGuess();
}

void ImplicitStepper::solveLinearFlow() {
  m_strand.dynamics().computeFlowMasses();

  computeFlowRHS();
  StrandDynamicTraits& dynamics = m_strand.dynamics();

  // precondition the velocity
  const int num_edges = m_strand.getNumEdges();
  const VecXx& mass = dynamics.getFlowMasses();
  for (int i = 0; i < num_edges; ++i) {
    if (mass(i) > 1e-20) {
      m_flow_newVelocities(i) = m_flow_rhs(i) / mass(i);
    } else {
      m_flow_newVelocities(i) = 0.0;
    }
  }

  computeFlowLHS();
  const ElasticStrandParameters& params = m_strand.getParameters();
  for (int i = 0; i < num_edges; ++i) {
    if (m_flow_lhs(i) > 1e-20) {
      m_flow_newVelocities(i) = m_flow_rhs(i) / m_flow_lhs(i);
    } else {
      m_flow_newVelocities(i) = 0.0;
    }

    // 3 CFL
    Scalar len = m_strand.getEdgeRestLength(i) * 3.0;
    if (fabs(m_flow_newVelocities(i)) > len / m_dt) {
      m_flow_newVelocities(i) = len / m_dt * sgn(m_flow_newVelocities(i));
    }
  }

  const int num_verts = m_strand.getNumVertices();
  const VecXx& flow_area = m_strand.getFutureState().getAreaDegreesOfFreedom();

  const Scalar N_times = m_strand.getParameters().getMaxFlowGradientRatio();

  // Gauss-Seidel iterations for smoothing
  for (int k = 0; k < 5; ++k) {
    for (int i = 1; i < num_verts - 1; ++i) {
      const Scalar a0 = flow_area(i - 1) + M_PI * m_strand.getRadiusA(i - 1) *
                                               m_strand.getRadiusB(i - 1);
      const Scalar a1 =
          flow_area(i) + M_PI * m_strand.getRadiusA(i) * m_strand.getRadiusB(i);
      const Scalar a2 = flow_area(i + 1) + M_PI * m_strand.getRadiusA(i + 1) *
                                               m_strand.getRadiusB(i + 1);

      if (a1 * 2.0 >= N_times * (a0 + a2)) {
        // remove coming flow
        const Scalar dv = m_flow_newVelocities(i - 1) - m_flow_newVelocities(i);
        const Scalar hv =
            (m_flow_newVelocities(i - 1) + m_flow_newVelocities(i)) * 0.5;
        if (dv > 0.) {
          m_flow_newVelocities(i) = m_flow_newVelocities(i - 1) = hv;
        }
      }
    }
  }

  // check_isnan("solve_linear_flow", m_flow_newVelocities);
}

void ImplicitStepper::computeGradAtVertex(const VecXx& edge_vars,
                                          VecXx& vertex_grads) {
  const int num_verts = m_strand.getNumVertices();
  const int num_edges = m_strand.getNumEdges();
  const ElasticStrandParameters& params = m_strand.getParameters();

  for (int i = 0; i < num_verts; ++i) {
    const int e_prev = i - 1;
    const int e_next = i;
    Scalar u_prev, u_next;

    if (e_prev >= 0) {
      u_prev = edge_vars(e_prev);
    } else {
      u_prev = edge_vars(e_next);
    }

    if (e_next < num_edges) {
      u_next = edge_vars(e_next);
    } else {
      u_next = edge_vars(e_prev);
    }

    vertex_grads(i) = (u_next - u_prev) / m_strand.getVoronoiLength(i);
  }

  // check_isnan("vertex_grads", vertex_grads);
}

void ImplicitStepper::accumulateReservoir() {
  const int idx_end = m_strand.getNumVertices() - 1;
  Vec2x& reservoir = m_strand.getReservoir();
  reservoir(0) =
      std::max(0.0, reservoir(0) + (m_strand.getFutureFlowDOFArea(0) +
                                    m_strand.getCurrentFlowDOFArea(0)) *
                                       m_flow_velocities(0) * m_dt * 0.5);
  reservoir(1) = std::max(
      0.0, reservoir(1) + (m_strand.getFutureFlowDOFArea(idx_end) +
                           m_strand.getCurrentFlowDOFArea(idx_end)) *
                              m_flow_velocities(idx_end - 1) * m_dt * 0.5);

  // check_isnan("reservoir_0", reservoir(0));
  // check_isnan("reservoir_1", reservoir(1));

  //        std::cout << "reservoir: " << reservoir << std::endl;
}

void ImplicitStepper::stepFlowAdvForce() {
  // ignore if no fluid is simulated
  if (!m_strand.getParent()->getFluidScriptingController()) return;

  // integrate force to velocity
  solveLinearFlow();
}

void ImplicitStepper::updateAdditionalInertia() {
  // NOTE: this must be called after backtracking!

  const int num_vtx = m_strand.getNumVertices();
  const int num_edges = m_strand.getNumEdges();

  std::vector<Scalar> summed_length(num_vtx);

  VecXx vert_vel(num_vtx);
  mapEdgeToVertices(m_flow_newVelocities, vert_vel);

  Scalar sum = 0.0;
  for (int i = 0; i < num_vtx; ++i) {
    summed_length[i] = sum;
    if (i < num_vtx - 1) sum += m_strand.getEdgeRestLength(i);
  }

  for (int i = 0; i < num_vtx; ++i) {
    Scalar pos = traceRK2Strand(summed_length[i], -m_dt, summed_length,
                                vert_vel, true, 0.0);
    if (i == num_vtx - 1) {
      const Scalar m = m_strand.getFutureSurfaceFlowMass(i);

      Vec4x u_s = interpolateStrandVelocity(
          pos, summed_length, m_savedVelocities, true, Vec4x::Zero());
      m_additional_inertia.segment<3>(i * 4) =
          (Vec3x(u_s(0), u_s(1), u_s(2)) -
           m_savedVelocities.segment<3>(i * 4)) *
          m;
    } else {
      const Scalar m = m_strand.getFutureSurfaceFlowMass(i);
      const Scalar m0 = m_strand.getCurrentSurfaceFlowMass(i);
      const Scalar I = m_strand.getEdgeFlowInertia(i);
      const Scalar I0 = m_strand.getCurrentEdgeFlowInertia(i);

      Vec4x u_s = interpolateStrandVelocity(
          pos, summed_length, m_savedVelocities, true, Vec4x::Zero());
      Mat4x mass_mat = Vec4x(m, m, m, I).asDiagonal();
      m_additional_inertia.segment<4>(i * 4) =
          mass_mat * (u_s - m_savedVelocities.segment<4>(i * 4));
    }
  }

  // check_isnan("updateAdditionalInertia", m_additional_inertia);
}

void ImplicitStepper::stepFlowData() {
  // ignore if no fluid is simulated
  if (!m_strand.getParent()->getFluidScriptingController()) return;
  // update flow data
  updateFlowData();

  accumulateReservoir();

  VecXx flow_area = m_strand.getFutureState().getAreaDegreesOfFreedom();

  // check_isnan("sfd_total_flow", m_total_flow);
  // check_isnan("area_step_flow_data_0", flow_area);

  Scalar new_total_flow = 0.0;
  for (int i = 0; i < flow_area.size() - 1; ++i) {
    new_total_flow +=
        (flow_area(i) + flow_area(i + 1)) * 0.5 * m_strand.getEdgeRestLength(i);
  }
  new_total_flow += m_strand.getReservoir()(0) + m_strand.getReservoir()(1);

  // check_isnan("sfd_total_flow", new_total_flow);

  if (new_total_flow > 1e-20) {
    for (int i = 0; i < flow_area.size(); ++i) {
      flow_area(i) *= m_total_flow / new_total_flow;
    }
  }

  // check_isnan("area_step_flow_data_1", flow_area);

  m_strand.getFutureState().setAreaDegreesOfFreedom(flow_area);

  m_strand.dynamics().invalidatePhysics();

  m_total_flow = new_total_flow;
}

void ImplicitStepper::updateFlowData() {
  const int num_verts = m_strand.getNumVertices();
  const ElasticStrandParameters& params = m_strand.getParameters();
  const auto controller = m_strand.getParent()->getFluidScriptingController();
  if (!controller) return;

  const int num_components = controller->getNumComponents();

  VecXx grad_u(num_verts);
  computeGradAtVertex(m_flow_newVelocities, grad_u);

  // check_isnan("area_update_grad_u", grad_u);

  VecXx flow_area = m_strand.getFutureState().getAreaDegreesOfFreedom();

  for (int i = 0; i < num_verts; ++i) {
    flow_area(i) *= exp(-grad_u(i) * m_dt);
    flow_area(i) = std::max(0.0, flow_area(i));

    const VecXx& color =
        m_flow_newComponents.segment(i * num_components, num_components);
    const Scalar mu = controller->getShearModulus(color);

    if (mu < 1e-20 || flow_area(i) < 1e-20) {
      m_flow_newStrain(i) = 0.0;
      continue;
    }

    m_flow_newStrain(i) =
        diffStrainElastic(m_flow_newStrain(i), 2.0 * grad_u(i), m_dt);

    // get parameter from color

    const Scalar mu2 = mu * mu;
    const Scalar yc = 0.8164965809 * controller->getYieldStress(color);
    const Scalar eps = controller->getCriterion();
    const int max_iters = controller->getMaxIters();
    const Scalar eta = controller->getViscosity(color);
    const Scalar n = controller->getFlowBehaviorIndex(color);
    const Scalar inv_n = 1.0 / n;
    const Scalar power_eta = pow(eta, inv_n);
    const Scalar coeff = 1.4142135624 * m_dt;

    // plasticity
    Scalar s_trial = 0.7071067812 * mu * fabs(m_flow_newStrain(i));
    Scalar phi_trial = s_trial - yc;
    if (phi_trial > 0.0) {
      Scalar sign_c = sgn(m_flow_newStrain(i));

      Scalar s_proj =
          bisection_root_finding(s_trial, yc, eps, max_iters, [&](Scalar s) {
            return power_eta * (s - s_trial) +
                   coeff * (s * sign_c + sqrt(s * s + 2.0 * mu2)) *
                       pow(s - yc, inv_n);
          });

      m_flow_newStrain(i) = 1.4142135624 * s_proj / mu * sign_c;
    }
  }

  m_strand.getFutureState().setAreaDegreesOfFreedom(flow_area);

  m_strand.dynamics().invalidatePhysics();

  // check_isnan("strain_update_flow_data", m_flow_newStrain);
  // check_isnan("area_update_flow_data", flow_area);

  //		std::cout << m_flow_newStrain << std::endl;
}

void ImplicitStepper::updateRHSwithImpulse(const VecXx& impulses) {
  m_rhs += impulses;

  StrandDynamicTraits& dynamics = m_strand.dynamics();
  dynamics.getScriptingController()->fixRHS(Lhs(), m_rhs, m_dt);
}

void ImplicitStepper::solveLinear() {
  StrandDynamicTraits& dynamics = m_strand.dynamics();

  computeRHS();
  computeLHS();

  Lhs().multiply(m_rhs, 1., m_newVelocities);
  dynamics.getScriptingController()->fixLHSAndRHS(Lhs(), m_rhs,
                                                  m_dt / m_fraction);

  m_linearSolver.store(Lhs());
}

void ImplicitStepper::setupFuturesFrames() {
  m_strand.getFutureState().m_referenceFrames1.set(
      m_strand.getCurrentState()
          .m_referenceFrames1
          .get());  // We need the old ones to compute the new ones
  m_strand.getFutureState().m_referenceFrames2.set(
      m_strand.getCurrentState()
          .m_referenceFrames2
          .get());  // We need the old ones to compute the new ones
  m_strand.getFutureState().m_referenceFrames1.getPreviousTangents() =
      m_strand.getCurrentState()
          .m_referenceFrames1
          .getPreviousTangents();  // Can we avoid the copis? // FIXME LATER
  m_strand.getFutureState().m_referenceTwists.set(
      m_strand.getCurrentState().m_referenceTwists.get());
}

void ImplicitStepper::computeRHS(bool dump_data, std::ostream& stream) {
  StrandDynamicTraits& dynamics = m_strand.dynamics();

  //         // check_isnan("rhs_0", m_rhs);

  m_rhs = m_velocities - m_newVelocities;

  //         // check_isnan("rhs_1", m_rhs);

  dynamics.multiplyByMassMatrix(m_rhs);

  //         // check_isnan("rhs_2", m_rhs);

  dynamics.computeFutureForces(!m_params.m_usePreFilterGeometry, true, false,
                               dump_data, stream);

  //        m_strand.getParameters().setKs( origKs );

  VecXx forces = m_strand.getFutureTotalForces();

  //         // check_isnan("rhs_forces", forces);

  m_rhs += forces * m_dt +
           (m_additional_impulse + m_additional_inertia) * m_fraction;

  if (dump_data) {
#pragma omp critical
    {
      stream << "[" << m_strand.getGlobalIndex() << "] RHS = " << m_rhs
             << std::endl;
    }
  }
}

Scalar ImplicitStepper::maxAdditionalImpulseNorm(int& idx) {
  Scalar maxlen = 0.;
  idx = -1;
  const int num_verts = m_strand.getNumVertices();
  for (int i = 0; i < num_verts; ++i) {
    const Scalar len = m_additional_impulse.segment<3>(i * 4).norm();
    //            maxlen = std::max(maxlen, len);
    if (len > maxlen) {
      maxlen = len;
      idx = i;
    }
  }
  return maxlen;
}

void ImplicitStepper::limitAdditionalImpulses(const Scalar& maxImpulses) {
  const int num_verts = m_strand.getNumVertices();
  for (int i = 0; i < num_verts; ++i) {
    const Scalar len = m_additional_impulse.segment<3>(i * 4).norm();
    if (len > maxImpulses) {
      m_additional_impulse.segment<3>(i * 4) =
          m_additional_impulse.segment<3>(i * 4).normalized() * maxImpulses;
    }
  }
}

VecXx& ImplicitStepper::impulse_rhs() {
  StrandDynamicTraits& dynamics = m_strand.dynamics();
  m_impulseRhs = m_newVelocities;
  dynamics.multiplyByMassMatrix(m_impulseRhs);

  if (m_params.m_useImpulseMethod) {
    dynamics.getScriptingController()->enforceVelocities(m_newVelocities, m_dt);
  } else {
    computeLHS();
    dynamics.getScriptingController()->fixLHSAndRHS(Lhs(), m_impulseRhs,
                                                    m_dt / m_fraction);
    m_linearSolver.store(Lhs());
  }

  return m_impulseRhs;
}

void ImplicitStepper::computeLHS(bool dump_data, std::ostream& stream) {
  StrandDynamicTraits& dynamics = m_strand.dynamics();

  dynamics.computeFutureJacobian(!m_params.m_usePreFilterGeometry, true, false,
                                 dump_data, stream);

  //        m_strand.getParameters().setKs( origKs );

  JacobianMatrixType& J = m_strand.getTotalJacobian();  // LHS = M - h^2 J
  J *= m_dt * m_dt;

  dynamics.addMassMatrixTo(J);

  if (dump_data) {
#pragma omp critical
    {
      stream << "[" << m_strand.getGlobalIndex() << "] LHS = " << J
             << std::endl;
    }
  }
}

void ImplicitStepper::computeFlowLHS() {
  StrandDynamicTraits& dynamics = m_strand.dynamics();
  m_flow_lhs = dynamics.getFlowMasses();

  if (!m_strand.getParent()) return;

  auto controller = m_strand.getParent()->getFluidScriptingController();
  if (!controller) return;

  const Scalar b = m_strand.getParameters().getSlipLength();

  const int num_components = controller->getNumComponents();
  const int num_edges = m_strand.getNumEdges();

  for (int i = 0; i < num_edges; ++i) {
    const Scalar avg_area = (m_strand.getFutureFlowDOFArea(i) +
                             m_strand.getFutureFlowDOFArea(i + 1)) *
                            0.5;
    const Scalar avg_radius_A =
        sqrt((m_strand.getRadiusA(i) * m_strand.getRadiusA(i) +
              m_strand.getRadiusA(i + 1) * m_strand.getRadiusA(i + 1)) *
             0.5);
    const Scalar avg_radius_B =
        sqrt((m_strand.getRadiusB(i) * m_strand.getRadiusB(i) +
              m_strand.getRadiusB(i + 1) * m_strand.getRadiusB(i + 1)) *
             0.5);

    const Scalar h = cyl_h_from_area(avg_radius_A, avg_radius_B, avg_area);
    const Scalar coeff = M_PI * (h + avg_radius_A + avg_radius_B) /
                         (b + h / 3.0) * m_strand.getEdgeRestLength(i);
    Scalar hu = h / std::max(fabs(m_flow_newVelocities(i)), 1e-7);

    VecXx color =
        m_flow_newComponents.segment(i * num_components, num_components) +
        m_flow_newComponents.segment((i + 1) * num_components, num_components);

    make_gibbs_simplex(color);

    const Scalar power = 1.0 - controller->getFlowBehaviorIndex(color);
    if (power < 0.0) hu = std::max(hu, 1e-10);

    const Scalar eta = controller->getViscosity(color);
    const Scalar tau_Y = controller->getYieldStress(color) * 0.8164965809;
    const Scalar tilde_eta = tau_Y * hu + eta * pow(hu, power);

    m_flow_lhs(i) += m_dt * coeff * tilde_eta;
  }

  // check_isnan("flow_lhs", m_flow_lhs);
}

void ImplicitStepper::markFlowStatus() {
  //		const int num_verts = m_strand.getNumVertices();
  //		const Scalar b = m_strand.getParameters().getSlipLength();
  //		m_flow_status.resize(num_verts);
  //
  //		for(int i = 0; i < num_verts; ++i) {
  //			const Scalar h = m_strand.getFutureFlowHeight(i);
  //			if(h < b) {
  //				m_flow_status[i] = 0U; // air
  //			} else {
  //				m_flow_status[i] = 3U; // liquid
  //			}
  //		}
  //
  //		for(int i = 0; i < num_verts; ++i) {
  //			if(m_flow_status[i] == 0U) {
  //				// ignore the case where both sides are liquid
  //				if(i > 0 && i < num_verts - 1 && m_flow_status[i - 1] == 3U
  //&& m_flow_status[i + 1] == 3U) { 					continue;
  //				}
  //
  //				if(i > 0 && m_flow_status[i - 1] == 3U) {
  //					m_flow_status[i] = 2U;
  //				} else if(i < num_verts - 1 && m_flow_status[i + 1] == 3U)
  //{ 					m_flow_status[i] = 1U;
  //				}
  //			}
  //		}
}

void ImplicitStepper::computeFlowRHS() {
  const auto controller = m_strand.getParent()->getFluidScriptingController();
  if (!controller) return;

  m_flow_rhs = m_flow_velocities;
  StrandDynamicTraits& dynamics = m_strand.dynamics();
  const int num_edges = m_strand.getNumEdges();
  const int num_verts = m_strand.getNumVertices();
  const int num_components = controller->getNumComponents();

  // gravity and fictitious force
  VecXx edge_accel(num_edges * 4);
  mapVertexToEdges<3, 4>(dynamics.getAccelerations(), edge_accel);

  for (int i = 0; i < num_edges; ++i) {
    m_flow_rhs(i) +=
        m_dt * (GravitationForce::getGravity() - edge_accel.segment<3>(i * 4))
                   .dot(m_strand.getFutureEdgeDirection(i));
  }

  dynamics.multiplyByFlowMatrix(m_flow_rhs);

  // shear stress
  for (int i = 0; i < num_edges; ++i) {
    const Scalar A0 = m_strand.getFutureFlowDOFArea(i);
    const Scalar A1 = m_strand.getFutureFlowDOFArea(i + 1);

    VecXx color =
        m_flow_newComponents.segment(i * num_components, num_components) +
        m_flow_newComponents.segment((i + 1) * num_components, num_components);

    make_gibbs_simplex(color);

    const Scalar shear_modulus = controller->getShearModulus(color);

    if (shear_modulus < 1e-20) continue;

    m_flow_rhs(i) += m_dt * shear_modulus * (A1 + A0) * 0.5 *
                     (m_flow_newStrain(i + 1) - m_flow_newStrain(i + 0));
  }

  // check_isnan("flow_rhs", m_flow_rhs);
}

const VecXx& ImplicitStepper::getDOFMasses() const {
  return m_strand.dynamics().getDOFMasses();
}

void ImplicitStepper::scale(const Scalar s) {
  if (s != 1.) {
    m_linearSolver.setScaling(s);
    m_rhs *= s;
  }
}

void ImplicitStepper::prepareDynamics() {
  // rest future state and intial guess
  m_strand.setFutureDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom());

  // update masses
  m_strand.dynamics().computeDOFMasses();

  if (m_params.m_useImpulseMethod) {
    // DK: for now init massMatrix here
    // (need to do this after computeViscousForceCoefficients since this where
    // mass is initialized for now -- awkward.)
    m_massMatrix.setZero();
    m_strand.dynamics().addMassMatrixTo(m_massMatrix);
    m_strand.dynamics().getScriptingController()->fixLHS(m_massMatrix);
    m_massMatrix_linearSolver.store(m_massMatrix);
  }

  if (m_params.m_usePreFilterGeometry) {
    filterGeometryLength(true);
    setupFuturesFrames();
  }

  StrandDynamicTraits& dynamics = m_strand.dynamics();

  m_newVelocities = dynamics.getDisplacements() / m_dt;
  m_strand.dynamics().getAccelerations() =
      (m_newVelocities - m_velocities) / m_dt;

  //		// check_isnan("newVels", m_newVelocities);

  dynamics.getDisplacements().setZero();
}

bool ImplicitStepper::prepareNewtonIteration() {
  StrandDynamicTraits& dynamics = m_strand.dynamics();

  dynamics.getDisplacements() = m_dt * m_newVelocities;
  const VecXx tentativeDofs =
      m_strand.getCurrentDegreesOfFreedom() + dynamics.getDisplacements();
  m_strand.setFutureDegreesOfFreedom(tentativeDofs);

  // check_isnan("pni_tentative", tentativeDofs);
  return false;
}

bool ImplicitStepper::performNewtonIteration() {
  if (m_newtonIter >= m_params.m_maxNewtonIterations) return true;

  StrandDynamicTraits& dynamics = m_strand.dynamics();
  const Scalar minAlpha = .01;  // Minimum step length

  m_prevRhs = m_rhs;
  computeRHS();

  // check_isnan("ni_rhs", rhs());

  if (m_newtonIter) {
    VecXx residual = m_rhs;

    dynamics.getScriptingController()->fixRHS(residual);
    const Scalar err = residual.squaredNorm() / residual.size();

    if (err < m_minErr) {
      m_minErr = err;

      // get it updated for possible later use
      m_bestLHS = Lhs();
      m_bestRhs = m_prevRhs;

      if (isSmall(err) || (m_newtonIter > 3 && m_minErr < 1e-6)) {
        m_rhs = m_prevRhs;
        return true;
      }
    }

    // Decrease or increase the step length based on current convergence
    if (err < m_prevErr) {
      m_alpha = std::min(1.0, 1.5 * m_alpha);
    } else {
      m_alpha = std::max(minAlpha, .5 * m_alpha);
    }

    m_prevErr = err;
  }

  //        std::cout << "NI[" << m_strand.getGlobalIndex() << "]: " <<
  //        m_newtonIter << ", " << m_prevErr << std::endl;

  computeLHS();

  // need to initialize them to avoid possible MKL error
  if (!m_newtonIter) {
    m_bestLHS = Lhs();
    m_bestRhs = m_rhs;
  }

  m_rhs = m_rhs * m_alpha;

  Lhs().multiply(m_rhs, 1., m_newVelocities);
  dynamics.getScriptingController()->fixLHSAndRHS(Lhs(), m_rhs,
                                                  m_dt / m_fraction);

  m_linearSolver.store(Lhs());
  m_notSPD = m_linearSolver.notSPD();
  m_linearSolver.solve(newVelocities(), rhs());

  // check_isnan("ni_velocity", m_newVelocities);

  m_strand.dynamics().getAccelerations() =
      (m_newVelocities - m_velocities) / m_dt;

  ++m_newtonIter;

  return false;
}

template <int nelem, int nstep>
void ImplicitStepper::mapVertexToEdges(const VecXx& vertex_vars,
                                       VecXx& edge_vars) {
  const int num_vtx = m_strand.getNumVertices();
  const int num_edges = m_strand.getNumEdges();

  if (num_edges == 0) return;

  for (int i = 0; i < num_edges; ++i) {
    edge_vars.template segment<nelem>(i * nstep) =
        (vertex_vars.template segment<nelem>(i * nstep) +
         vertex_vars.template segment<nelem>((i + 1) * nstep)) *
        0.5;
  }
}

template <int nelem, int nstep>
void ImplicitStepper::mapEdgeToVertices(const VecXx& edge_vars,
                                        VecXx& vertex_vars) {
  const int num_vtx = m_strand.getNumVertices();
  const int num_edges = m_strand.getNumEdges();
  if (num_edges == 0) return;

  for (int i = 0; i < num_vtx; ++i) {
    if (i == 0) {
      vertex_vars.template segment<nelem>(0) =
          edge_vars.template segment<nelem>(0);
    } else if (i == num_vtx - 1) {
      vertex_vars.template segment<nelem>(i * nstep) =
          edge_vars.template segment<nelem>((num_edges - 1) * nstep);
    } else {
      const int e_prev = i - 1;
      const int e_next = i;

      vertex_vars.template segment<nelem>(i * nstep) =
          (edge_vars.template segment<nelem>(e_prev * nstep) +
           edge_vars.template segment<nelem>(e_next * nstep)) *
          0.5;
    }
  }
}

Vec4x ImplicitStepper::interpolateStrandVelocity(
    const Scalar& pos, const std::vector<Scalar>& summed_length,
    const VecXx& source, bool clamped, const Vec4x& default_val) {
  Scalar pos_strand = find_in_ordered(summed_length, pos);
  int np = source.size() / 4;
  int ip = (int)floor(pos_strand);
  if (clamped && (ip < 0 || ip >= np - 1)) {
    return default_val;
  }

  ip = clamp(ip, 0, np - 2);
  Scalar frac = clamp(pos_strand - (Scalar)ip, 0.0, 1.0);
  Vec4x s0 = source.segment<4>(ip * 4);
  Vec4x s1;
  if (ip + 1 == np - 1) {
    s1 = Vec4x(source((ip + 1) * 4 + 0), source((ip + 1) * 4 + 1),
               source((ip + 1) * 4 + 2), source(ip * 4 + 3));
  } else {
    s1 = source.segment<4>((ip + 1) * 4);
  }

  return lerp(s0, s1, frac);
}

Scalar ImplicitStepper::interpolateStrand(
    const Scalar& pos, const std::vector<Scalar>& summed_length,
    const VecXx& source, bool clamped, const Scalar& default_val) {
  Scalar pos_strand = find_in_ordered(summed_length, pos);
  int np = source.size();
  int ip = (int)floor(pos_strand);
  if (clamped && (ip < 0 || ip >= np - 1)) {
    return default_val;
  }

  ip = clamp(ip, 0, np - 2);
  Scalar frac = clamp(pos_strand - (Scalar)ip, 0.0, 1.0);
  return lerp(source[ip], source[ip + 1], frac);
}

VecXx ImplicitStepper::interpolateStrand(
    const Scalar& pos, const std::vector<Scalar>& summed_length,
    const VecXx& source, bool clamped, const VecXx& default_val) {
  const int num_components = default_val.size();
  Scalar pos_strand = find_in_ordered(summed_length, pos);
  int np = source.size() / num_components;
  int ip = (int)floor(pos_strand);
  if (clamped && (ip < 0 || ip >= np - 1)) {
    return default_val;
  }

  ip = clamp(ip, 0, np - 2);
  Scalar frac = clamp(pos_strand - (Scalar)ip, 0.0, 1.0);
  VecXx ret =
      source.segment(ip * num_components, num_components) * (1.0 - frac) +
      source.segment((ip + 1) * num_components, num_components) * frac;
  return ret;
}

Scalar ImplicitStepper::traceRK2Strand(const int ipos, const Scalar& dt,
                                       const std::vector<Scalar>& summed_length,
                                       const VecXx& u_vert, bool clamped,
                                       const Scalar& default_val) {
  Scalar p = summed_length[ipos];
  Scalar vel = u_vert[ipos];
  vel = interpolateStrand(p + vel * dt * 0.5, summed_length, u_vert, clamped,
                          vel);
  p += dt * vel;
  return p;
}

Scalar ImplicitStepper::traceRK2Strand(const Scalar& pos, const Scalar& dt,
                                       const std::vector<Scalar>& summed_length,
                                       const VecXx& u_vert, bool clamped,
                                       const Scalar& default_val) {
  Scalar p = pos;
  Scalar vel =
      interpolateStrand(p, summed_length, u_vert, clamped, default_val);
  vel = interpolateStrand(p + vel * dt * 0.5, summed_length, u_vert, clamped,
                          vel);
  p += dt * vel;
  return p;
}

void ImplicitStepper::backtrackFlowData() {
  const int num_vtx = m_strand.getNumVertices();
  const int num_edges = m_strand.getNumEdges();
  StrandState& state = m_strand.getFutureState();

  std::vector<Scalar> summed_length(num_vtx);

  VecXx vert_vel(num_vtx);
  mapEdgeToVertices(m_flow_newVelocities, vert_vel);

  Scalar sum = 0.0;
  for (int i = 0; i < num_vtx; ++i) {
    summed_length[i] = sum;
    if (i < num_vtx - 1) sum += m_strand.getEdgeRestLength(i);
  }

  for (int i = 0; i < num_edges; ++i) {
    Scalar pos = (summed_length[i] + summed_length[i + 1]) * 0.5;
    pos = traceRK2Strand(pos, -m_dt, summed_length, vert_vel, true, 0.0);
    m_flow_newVelocities(i) =
        interpolateStrand(pos, summed_length, vert_vel, true, 0.0);
  }

  VecXx newFlowArea(num_vtx);
  const VecXx& curFlowArea =
      m_strand.getCurrentState().getAreaDegreesOfFreedom();

  // record total volume
  m_total_flow = 0.0;
  for (int i = 0; i < curFlowArea.size() - 1; ++i) {
    m_total_flow += (curFlowArea(i) + curFlowArea(i + 1)) * 0.5 *
                    m_strand.getEdgeRestLength(i);
  }
  m_total_flow += m_strand.getReservoir()(0) + m_strand.getReservoir()(1);

  // check_isnan("backtrack_total_flow", m_total_flow);

  const int num_components = m_flow_components.size() / num_vtx;
  VecXx default_comps(num_components);
  default_comps.setZero();
  default_comps(0) = 1.0;

  // first compute A * c
  VecXx ac(m_flow_components.size());
  for (int i = 0; i < num_vtx; ++i) {
    ac.segment(i * num_components, num_components) =
        m_flow_components.segment(i * num_components, num_components) *
        curFlowArea(i);
  }

  //        std::cout << ac << std::endl;

  for (int i = 0; i < num_vtx; ++i) {
    Scalar pos = traceRK2Strand(i, -m_dt, summed_length, vert_vel, true, 0.0);
    newFlowArea(i) =
        interpolateStrand(pos, summed_length, curFlowArea, true, 0.0);
    m_flow_newStrain(i) =
        interpolateStrand(pos, summed_length, m_flow_strain, true, 0.0);
    m_flow_newComponents.segment(i * num_components, num_components) =
        interpolateStrand(pos, summed_length, ac, false, default_comps);
    Scalar sum =
        m_flow_newComponents.segment(i * num_components, num_components).sum();
    if (sum < 1e-20) {
      m_flow_newComponents.segment(i * num_components, num_components) =
          default_comps;
    } else {
      m_flow_newComponents.segment(i * num_components, num_components) /= sum;
    }
  }

  //        std::cout << m_flow_newComponents << std::endl;

  state.setAreaDegreesOfFreedom(newFlowArea);

  m_strand.dynamics().invalidatePhysics();

  // check_isnan("backtrack_flow_data_v", m_flow_newVelocities);
  // check_isnan("backtrack_flow_data_a", newFlowArea);
  // check_isnan("backtrack_flow_data_c", m_flow_newStrain);
}

void ImplicitStepper::prepareSolveNonlinear() {
  m_minErr = 1e+99;
  m_prevErr = 1e+99;
  m_newtonIter = 0;

  m_alpha = 1.0;

  m_strand.requireExactJacobian(m_params.m_useExactJacobian);
  m_strand.projectJacobian(m_params.m_useProjectedJacobian);
}

bool ImplicitStepper::postSolveNonlinear() {
  if (m_newtonIter == m_params.m_maxNewtonIterations) {
    // Failed. We need to rollback with lower stretch multiplier
    //            m_stretchMultiplier *= 0.5;
    //
    //            StrandDynamicTraits& dynamics = m_strand.dynamics() ;
    //
    //            dynamics.getDisplacements() = m_dt * m_velocities ;
    //
    //            return false;

    // std::cout << m_newtonIter << std::endl;
    m_rhs = m_bestRhs;
    Lhs() = m_bestLHS;

    m_linearSolver.store(Lhs());
    m_linearSolver.solve(newVelocities(), rhs());
    m_notSPD = m_linearSolver.notSPD();

    m_strand.dynamics().getAccelerations() =
        (m_newVelocities - m_velocities) / m_dt;

    // Dump Data
    //            prepareNewtonIteration();
    //
    //            const std::string& output_dir =
    //            m_strand.getParent()->getOutputDirectory(); Scalar time =
    //            m_strand.getParent()->getTime(); int global_idx =
    //            m_strand.getGlobalIndex();
    //
    //            std::ostringstream oss;
    //            oss << output_dir << "/" << "dump_" << global_idx << "_" <<
    //            time << ".log";
    //
    //            std::ofstream ofs(oss.str().c_str());
    //
    //            computeLHS(true, ofs);
    //            computeRHS(true, ofs);
    //
    //            ofs.close();
  }

  m_usedNonlinearSolver = true;

  return true;

  // check_isnan("psn_velocity", m_newVelocities);

  //        std::cout << "Nonlinear solve for strand " <<
  //        m_strand.getGlobalIndex() << " iter " << m_newtonIter << " err " <<
  //        m_minErr << " ; SPD: " << (int) !m_notSPD << std::endl ;
}

// Updates the current Lhs and rhs based of the m_newVelocities guess
bool ImplicitStepper::updateLinearSystem(const VecXx solverForces)
// DK: this is where the lineic stretch gets checked and a possible update is
// applied
{
  StrandDynamicTraits& dynamics = m_strand.dynamics();

  m_strand.requireExactJacobian(false);

  const VecXx tentativeDofs =
      m_strand.getCurrentDegreesOfFreedom() + m_dt * m_newVelocities;
  m_strand.setFutureDegreesOfFreedom(tentativeDofs);

  dynamics.getDisplacements() = m_dt * m_newVelocities;

  const Scalar stretchE = getLineicStretch();

  bool needUpdate = stretchE > m_stretchingFailureThreshold;

  ContactStream(g_log, "GS")
      << "Strand " << m_strand.getGlobalIndex() << " stretch is : " << stretchE
      << " ( with max stretch : " << m_stretchingFailureThreshold << " )";

  if (!needUpdate) {
    computeRHS();
    dynamics.getScriptingController()->fixRHS(m_rhs);

    m_usedNonlinearSolver = true;
    const Scalar residual = (m_rhs + solverForces).squaredNorm() / m_rhs.rows();

    needUpdate = residual > m_costretchResidualFailureThreshold;

    ContactStream(g_log, "GS")
        << "Strand " << m_strand.getGlobalIndex() << " residual is " << residual
        << " ( with max residual : " << m_costretchResidualFailureThreshold
        << " )";
  }

  ContactStream(g_log, "GS") << "Strand " << m_strand.getGlobalIndex()
                             << " needUpdate = " << needUpdate;

  if (needUpdate)  // DK: if residual or stretchE too big.
  {
    solveLinear();
    // m_strand.swapStates();
    // dynamics.flashForRendering() ;
    // m_strand.swapStates();
  }

  return needUpdate;  // DK: now returns true if needs update
}

Scalar ImplicitStepper::getLineicStretch() {
  Scalar stretchE = 0;
  ForceAccumulator<StretchingForce<NonViscous> >::accumulateFuture(stretchE,
                                                                   m_strand);
  stretchE /= m_strand.getParameters().getKs() * m_strand.getTotalRestLength();

  return stretchE;
}

// Updates the strand using m_newVelocities and check for stretching. If
// stretching, calls the appropriate fail-safe or returns false DK: here's where
// most of the failsafes kick in -- should modify here to best test algorithm
//     note that this is called by
//     StrandImplicitManager::postProcessFrictionProblem() where the bool
//     returned here is ignored but m_lastStepWasRejected is set which is
//     checked in StrandImplicitManager::solveCollidingGroup() to see if
//     failsafe is needed
bool ImplicitStepper::update(bool afterConstraints) {
  VecXx displacements = m_newVelocities * m_dt;

  // check_isnan("displacements in update", displacements);

  m_strand.dynamics().getScriptingController()->enforceDisplacements(
      displacements, m_fraction);

  m_strand.setFutureDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom() +
                                     displacements);

  bool accept = true;

  // FTL projection here as strong failsafes, we project all strands that
  // reaches maximal # iters
  if (m_params.m_useLengthProjection) {
    m_strand.dynamics().getScriptingController()->enforceDisplacements(
        displacements);

    m_newVelocities = displacements / m_dt;

    m_strand.setFutureDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom() +
                                       displacements);

    filterGeometryLength(false);
  } else if (!m_params.m_usePreFilterGeometry) {
    // DK: failsafes here:
    const Scalar stretchE = getLineicStretch();
    if (stretchE > m_stretchingFailureThreshold) {
      //  Don't bother with very small strands, they would create numerical
      //  problems anyway
      if (m_strand.getNumVertices() == 3 &&
          m_strand.getTotalRestLength() < .1) {
        WarningStream(g_log, "") << "Strand " << m_strand.getGlobalIndex()
                                 << " did bad stuff and is small (len="
                                 << m_strand.getTotalRestLength()
                                 << "), reverting to unsimulated pos";

        VecXx futureDofs(m_strand.getCurrentDegreesOfFreedom());
        //                    m_strand.getScriptingController()->computeRigidBodyMotion(
        //                                futureDofs,
        //                                m_strand.getCurrentDegreesOfFreedom()
        //                                );
        m_strand.setFutureDegreesOfFreedom(futureDofs);
        m_strand.filterFutureGeometryByRestLength(1.e-6);

      } else if (afterConstraints) {
        accept = false;
      }
    }
  }
  m_lastStepWasRejected = !accept;

  m_strand.dynamics()
      .acceptGuess();  // DK: here's where the solution gets transfered from
                       // futurestate to currentstate (note *not* in finalize())
                       // and then current state is possibly accepted as "good"

  m_strand.getFutureState().freeCachedQuantities();

  m_savedVelocities = m_newVelocities;

  return accept;
}

// DK: relax thetas in finalize is yet another failsafe here...
void ImplicitStepper::finalize() {
  if (!m_params.m_simulationManager_limitedMemory) {
    m_strand.dynamics().flashForRendering();
  }
  m_strand.dynamics().nanFailSafe();

  updateRodAccelerationAndVelocity();

  m_flow_velocities = m_flow_newVelocities;

  m_flow_strain = m_flow_newStrain;

  m_flow_components = m_flow_newComponents;
  // update multipliers
  //    m_strand.accumulateMultipliers(true, false);
}

void ImplicitStepper::updateRodAcceleration() {
  m_strand.dynamics().getAccelerations() =
      (m_newVelocities - m_velocities) / m_dt;
}

void ImplicitStepper::updateRodAccelerationAndVelocity() {
  // Finite difference computation of velocities and acceleration
  newVelocities() = m_strand.dynamics().getDisplacements() / m_dt;
  m_strand.dynamics().getAccelerations() =
      (m_newVelocities - m_velocities) / m_dt;

  // Here's where the old velocities are updated.
  m_velocities = m_newVelocities;
}

void ImplicitStepper::rewind() { m_strand.swapStates(); }

void ImplicitStepper::filterGeometryLength(bool preStep) {
  //#pragma omp critical
  //        {
  //            std::cout << "[filter strand " << m_strand.getGlobalIndex() <<
  //            "]" << std::endl;
  //        }
  //
  VecXx futureDofs = m_strand.getFutureDegreesOfFreedom();

  m_strand.filterFutureGeometryByRestLength(0.0, false);

  m_projectionDisplacements = m_strand.getFutureDegreesOfFreedom() - futureDofs;

  m_velocities = (m_strand.getFutureDegreesOfFreedom() -
                  m_strand.getCurrentDegreesOfFreedom()) /
                 m_dt;

  // correct velocities
  const int num_verts = m_strand.getNumVertices();
  for (int i = 0; i < num_verts; ++i) {
    if (m_strand.isVertexFreezed(i) || m_strand.isVertexGoaled(i)) continue;

    if (i < num_verts - 1)
      m_newVelocities.segment<3>(i * 4) +=
          (m_projectionDisplacements.segment<3>(i * 4) -
           m_projectionDisplacements.segment<3>((i + 1) * 4)) /
          m_dt;
    else
      m_newVelocities.segment<3>(i * 4) +=
          m_projectionDisplacements.segment<3>(i * 4) / m_dt;
  }

  m_strand.dynamics().getDisplacements() = m_newVelocities * m_dt;
}

const JacobianMatrixType& ImplicitStepper::Lhs() const {
  return m_strand.getTotalJacobian();
}

JacobianMatrixType& ImplicitStepper::Lhs() {
  return m_strand.getTotalJacobian();
}

}  // namespace strandsim
