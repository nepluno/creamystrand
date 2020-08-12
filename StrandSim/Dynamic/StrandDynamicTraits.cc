/**
 * \copyright 2014 Danny Kaufman
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "StrandDynamicTraits.hh"

#include "../Core/ElasticStrand.hh"
#include "../Core/ElasticStrandUtils.hh"
#include "../Forces/AirDragForce.hh"
#include "../Forces/BendingForce.hh"
#include "../Forces/FixingForce.hh"
#include "../Forces/FluidDragForce.hh"
#include "../Forces/FluidPressureForce.hh"
#include "../Forces/ForceAccumulator.hh"
#include "../Forces/GravitationForce.hh"
#include "../Forces/StretchingForce.hh"
#include "../Forces/TwistingForce.hh"
#include "../Forces/ViscousOrNotViscous.hh"
#include "../Render/StrandRenderer.hh"
#include "../Utils/MathUtilities.hh"
#include "../Utils/TextLog.hh"
#include "DOFScriptingController.hh"
#include "FluidScriptingController.hh"
#include "StrandImplicitManager.hh"

namespace strandsim {

bool StrandDynamicTraits::s_freeMemoryAfterComputeJacobian = true;

StrandDynamicTraits::StrandDynamicTraits(ElasticStrand& strand)
    : m_strand(strand),
      m_scriptingController(NULL),      //
      m_futureJacobianUpToDate(false),  //
      m_futureForcesUpToDate(false),    //
      m_DOFmassesUpToDate(false),
      m_flowMassesUpToDate(false),
      m_renderer(NULL) {}

StrandDynamicTraits::~StrandDynamicTraits() {}

void StrandDynamicTraits::resizeSelf() {
  const unsigned ndofs = m_strand.getCurrentDegreesOfFreedom().rows();

  m_DOFmasses.resize(ndofs);
  m_flowMasses.resize(m_strand.getNumEdges());

  m_displacements.resize(ndofs);
  m_displacements.setZero();
  m_accelerations.resize(ndofs);
  m_accelerations.setZero();
  m_immuneEdges.resize(m_strand.getNumEdges(), false);
  m_outsideVertices.resize(m_strand.getNumVertices(), true);
  m_nearVertices.resize(m_strand.getNumVertices(), false);
  m_insideEdges.resize(0);
  m_interfacingEdges.resize(m_strand.getNumEdges(), false);
  m_passingVertices.resize(m_strand.getNumVertices(), false);
}

void StrandDynamicTraits::computeViscousForceCoefficients(Scalar dt) {
  m_strand.getParameters().computeViscousForceCoefficients(dt);
}

void StrandDynamicTraits::computeDOFMasses() {
  if (m_DOFmassesUpToDate) return;

  for (IndexType vtx = 0; vtx < m_strand.m_numVertices; ++vtx) {
    m_DOFmasses[4 * vtx + 0] = m_DOFmasses[4 * vtx + 1] =
        m_DOFmasses[4 * vtx + 2] = m_strand.m_vertexMasses[vtx] +
                                   m_strand.getFutureSurfaceFlowMass(vtx);

    if (vtx < m_strand.m_numEdges)
      m_DOFmasses[4 * vtx + 3] = m_strand.getEdgeInertia(vtx);
  }

  m_DOFmassesUpToDate = true;
}

void StrandDynamicTraits::computeFlowMasses() {
  if (m_flowMassesUpToDate) return;

  for (IndexType edge = 0; edge < m_strand.m_numEdges; ++edge) {
    m_flowMasses[edge] = m_strand.getFutureSurfaceFlowMassAtEdge(edge);
  }

  m_flowMassesUpToDate = true;
}

const VecXx& StrandDynamicTraits::getFlowMasses() const { return m_flowMasses; }

void StrandDynamicTraits::multiplyByFlowMatrix(VecXx& F) const {
  F.array() *= m_flowMasses.array();
}
////////////////////////////////////////////////////////////////////////////////
// Dynamic methods, using viscous forces
////////////////////////////////////////////////////////////////////////////////

void StrandDynamicTraits::computeFutureJacobian(bool withStretch,
                                                bool withViscous,
                                                bool butOnlyForBendingModes,
                                                bool dump_data,
                                                std::ostream& dump_stream) {
  if (m_futureJacobianUpToDate) {
    return;
  }

  StrandState& futureState = *m_strand.m_futureState;

  JacobianMatrixType& futureJ = *(futureState.m_totalJacobian);
  futureJ.setZero();

  m_strand.accumulateJ<FixingForce<NonViscous> >(futureState);

  if (dump_data) {
#pragma omp critical
    {
      dump_stream << "[" << m_strand.getGlobalIndex() << "] FJ0 = " << futureJ
                  << std::endl;
    }
  }

  if (withStretch) {
    m_strand.accumulateJ<StretchingForce<NonViscous> >(futureState);
    if (dump_data) {
#pragma omp critical
      {
        dump_stream << "[" << m_strand.getGlobalIndex() << "] FJ1 = " << futureJ
                    << std::endl;
      }
    }
  }

  m_strand.accumulateJ<TwistingForce<NonViscous> >(futureState);

  if (dump_data) {
#pragma omp critical
    {
      dump_stream << "[" << m_strand.getGlobalIndex() << "] FJ2 = " << futureJ
                  << std::endl;
    }
  }

  m_strand.accumulateJ<BendingForce<NonViscous> >(futureState);

  if (dump_data) {
#pragma omp critical
    {
      dump_stream << "[" << m_strand.getGlobalIndex() << "] FJ3 = " << futureJ
                  << std::endl;
    }
  }

  if (withViscous) {
    if (!butOnlyForBendingModes) {
      if (withStretch) {
        m_strand.accumulateJ<StretchingForce<Viscous> >(futureState);

        if (dump_data) {
#pragma omp critical
          {
            dump_stream << "[" << m_strand.getGlobalIndex()
                        << "] FJ5 = " << futureJ << std::endl;
          }
        }
      }
    }
    m_strand.accumulateJ<TwistingForce<Viscous> >(futureState);
    if (dump_data) {
#pragma omp critical
      {
        dump_stream << "[" << m_strand.getGlobalIndex() << "] FJ6 = " << futureJ
                    << std::endl;
      }
    }

    m_strand.accumulateJ<BendingForce<Viscous> >(futureState);

    if (dump_data) {
#pragma omp critical
      {
        dump_stream << "[" << m_strand.getGlobalIndex() << "] FJ7 = " << futureJ
                    << std::endl;
      }
    }

    m_strand.accumulateJ<AirDragForce>(futureState);
    if (dump_data) {
#pragma omp critical
      {
        dump_stream << "[" << m_strand.getGlobalIndex() << "] FJ8 = " << futureJ
                    << std::endl;
      }
    }
    m_strand.accumulateJ<FluidDragForce>(futureState);
    if (dump_data) {
#pragma omp critical
      {
        dump_stream << "[" << m_strand.getGlobalIndex() << "] FJ9 = " << futureJ
                    << std::endl;
      }
    }
  }

  m_strand.accumulateJ<FluidPressureForce>(futureState);
  if (dump_data) {
#pragma omp critical
    {
      dump_stream << "[" << m_strand.getGlobalIndex() << "] FJ10 = " << futureJ
                  << std::endl;
    }
  }

  m_strand.accumulateRuntimeForces(RuntimeForceBase::J, futureState);

  futureJ *= -1.0;  // To match BASim's sign conventions

  m_futureJacobianUpToDate = true;

  if (dump_data) {
#pragma omp critical
    {
      dump_stream << "[" << m_strand.getGlobalIndex() << "] FJ13 = " << futureJ
                  << std::endl;
    }
  }

  // Free some memory
  if (s_freeMemoryAfterComputeJacobian) {
    futureState.m_hessTwists.free();
    futureState.m_hessKappas.free();
  }
}

void StrandDynamicTraits::computeLHS(Scalar dt, bool withStretch,
                                     bool withViscous, bool dump_data,
                                     std::ostream& dump_stream) {
  computeFutureJacobian(withStretch, withViscous);
  JacobianMatrixType& LHS = m_strand.getTotalJacobian();
  LHS *= dt * dt;
  addMassMatrixTo(LHS);
  getScriptingController()->fixLHS(LHS);  // Enforce scripted vertices

  if (dump_data) {
#pragma omp critical
    {
      dump_stream << "[" << m_strand.getGlobalIndex() << "] LHS = " << LHS
                  << std::endl;
    }
  }
}

void StrandDynamicTraits::computeFutureForces(bool withStretch,
                                              bool withViscous,
                                              bool butOnlyForBendingModes,
                                              bool dump_data,
                                              std::ostream& dump_stream) {
  if (m_futureForcesUpToDate) {
    return;
  }

  StrandState& futureState = *m_strand.m_futureState;

  futureState.m_totalEnergy = 0.0;  // NB energy is not going to be used
  VecXx& futureF = futureState.m_totalForce;
  futureF.setZero();

  m_strand.accumulateF<FixingForce<NonViscous> >(futureState);

  if (dump_data) {
#pragma omp critical
    {
      dump_stream << "[" << m_strand.getGlobalIndex()
                  << "] FF0 = " << futureF.transpose() << std::endl;
    }
  }

  if (withStretch) {
    m_strand.accumulateF<StretchingForce<NonViscous> >(futureState);
    // check_isnan("force_0", futureF);

    if (dump_data) {
#pragma omp critical
      {
        dump_stream << "[" << m_strand.getGlobalIndex()
                    << "] FF1 = " << futureF.transpose() << std::endl;
      }
    }
  }

  m_strand.accumulateF<TwistingForce<NonViscous> >(futureState);
  // check_isnan("force_1", futureF);
  if (dump_data) {
#pragma omp critical
    {
      dump_stream << "[" << m_strand.getGlobalIndex()
                  << "] FF2 = " << futureF.transpose() << std::endl;
    }
  }

  m_strand.accumulateF<BendingForce<NonViscous> >(futureState);
  // check_isnan("force_2", futureF);
  if (dump_data) {
#pragma omp critical
    {
      dump_stream << "[" << m_strand.getGlobalIndex()
                  << "] FF3 = " << futureF.transpose() << std::endl;
    }
  }
  if (withViscous) {
    if (!butOnlyForBendingModes) {
      if (withStretch) {
        m_strand.accumulateF<StretchingForce<Viscous> >(futureState);
        if (dump_data) {
#pragma omp critical
          {
            dump_stream << "[" << m_strand.getGlobalIndex()
                        << "] FF5 = " << futureF.transpose() << std::endl;
          }
        }
      }
      // check_isnan("force_3", futureF);
    }
    m_strand.accumulateF<TwistingForce<Viscous> >(futureState);
    if (dump_data) {
#pragma omp critical
      {
        dump_stream << "[" << m_strand.getGlobalIndex()
                    << "] FF6 = " << futureF.transpose() << std::endl;
      }
    }
    // check_isnan("force_4", futureF);

    m_strand.accumulateF<BendingForce<Viscous> >(futureState);
    if (dump_data) {
#pragma omp critical
      {
        dump_stream << "[" << m_strand.getGlobalIndex()
                    << "] FF7 = " << futureF.transpose() << std::endl;
      }
    }
    // check_isnan("force_5", futureF);

    m_strand.accumulateF<AirDragForce>(futureState);
    if (dump_data) {
#pragma omp critical
      {
        dump_stream << "[" << m_strand.getGlobalIndex()
                    << "] FF8 = " << futureF.transpose() << std::endl;
      }
    }
    // check_isnan("force_6", futureF);

    m_strand.accumulateF<FluidDragForce>(futureState);
    if (dump_data) {
#pragma omp critical
      {
        dump_stream << "[" << m_strand.getGlobalIndex()
                    << "] FF9 = " << futureF.transpose() << std::endl;
      }
    }
    // check_isnan("force_7", futureF);
  }
  // check_isnan("force_8", futureF);

  m_strand.accumulateF<GravitationForce>(futureState);
  if (dump_data) {
#pragma omp critical
    {
      dump_stream << "[" << m_strand.getGlobalIndex()
                  << "] FF11 = " << futureF.transpose() << std::endl;
    }
  }
  // check_isnan("force_9", futureF);

  m_strand.accumulateF<FluidPressureForce>(futureState);
  if (dump_data) {
#pragma omp critical
    {
      dump_stream << "[" << m_strand.getGlobalIndex()
                  << "] FF12 = " << futureF.transpose() << std::endl;
    }
  }
  // check_isnan("force_10", futureF);

  m_strand.accumulateRuntimeForces(RuntimeForceBase::EF, futureState);
  // check_isnan("force_12", futureF);

  m_futureForcesUpToDate = true;
}

void StrandDynamicTraits::computeFutureConservativeEnergy(bool withStretch) {
  StrandState& futureState = *m_strand.m_futureState;

  if (withStretch)
    m_strand.accumulateE<StretchingForce<NonViscous> >(futureState);

  m_strand.accumulateE<TwistingForce<NonViscous> >(futureState);
  m_strand.accumulateE<BendingForce<NonViscous> >(futureState);
  m_strand.accumulateE<GravitationForce>(futureState);

  m_strand.accumulateRuntimeForces(RuntimeForceBase::E, futureState);
}

void StrandDynamicTraits::addMassMatrixTo(JacobianMatrixType& J,
                                          Scalar multiplier) const {
  for (int i = 0; i < m_DOFmasses.size(); i++) {
    J(i, i) += m_DOFmasses[i] * multiplier;
  }
}

const VecXx& StrandDynamicTraits::getDOFMasses() const { return m_DOFmasses; }

void StrandDynamicTraits::multiplyByMassMatrix(VecXx& tmp,
                                               Scalar multiplier) const {
  tmp.array() *= (m_DOFmasses * multiplier).array();
}

void StrandDynamicTraits::acceptGuess() {
  // If we are here, that's because we exited the Newton loop before futureState
  // was touched, so if we swap currentState will contain the accepted guess and
  // the up-to-date forces and Jacobian.

  m_displacements = m_strand.getFutureDegreesOfFreedom() -
                    m_strand.getCurrentDegreesOfFreedom();

  m_strand.swapStates();
}

void StrandDynamicTraits::acceptDisplacements() {
  m_strand.setFutureDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom() +
                                     m_displacements);

  m_strand.swapStates();
}

int StrandDynamicTraits::flashForRendering() {
  m_renderer->m_debugVertices.push_back(std::vector<float>());

  std::vector<float>& currentVertices = m_renderer->m_debugVertices.back();
  currentVertices.resize(m_strand.getNumVertices() * 3);
  Eigen::Map<Eigen::VectorXf> vert(&currentVertices[0], currentVertices.size());
  m_strand.getCurrentVertices(vert);

  return m_renderer->m_debugVertices.size();
}

bool StrandDynamicTraits::isImmune(int edge) const {
  return m_immuneEdges[edge];
}

void StrandDynamicTraits::addHitPoint(const Vec3x& point) {
  m_renderer->m_hitPoints.push_back(point.cast<float>());
}

void StrandDynamicTraits::clearDebugDrawing() {
  m_renderer->m_debugVertices.clear();
  m_renderer->m_hitPoints.clear();
}

void StrandDynamicTraits::nanFailSafe() {
  if (containsNans(m_strand.getCurrentDegreesOfFreedom())) {
    ErrorStream(g_log, "")
        << "Elastic strand " << m_strand.m_globalIndex
        << " was contaminated by NaNs: reverting to rigid motion";
    m_strand.swapStates();

    const VecXx& currentDOFs = m_strand.getCurrentDegreesOfFreedom();
    VecXx futureDOFs;
    m_scriptingController->computeRigidBodyMotion(futureDOFs, currentDOFs);

    if (containsNans(futureDOFs)) {
      ErrorStream(g_log, "") << "Damn it. ";
      exit(-1);
    }

    m_strand.setFutureDegreesOfFreedom(futureDOFs);

    acceptGuess();  // Re-compute displacements as well, as they were problably
                    // NaNised
  }
}

void StrandDynamicTraits::checkPassingVertices() {
  auto controller = m_strand.getParent()->getFluidScriptingController();
  if (!controller) return;

  const int num_verts = m_strand.getNumVertices();
  const Scalar dx = controller->getDx();

  m_passingVertices.assign(num_verts, false);

  // check and buffer all vertices phi
  std::vector<Scalar> distances(num_verts);
  for (int i = 0; i < num_verts; ++i) {
    const Vec3x pos = m_strand.getVertex(i);
    distances[i] = controller->get_dense_liquid_signed_distance(pos);
  }

  for (int i : m_promisingVertices) {
    // vertex passed through interface during last time step
    m_passingVertices[i] = distances[i] >= -0.5 * dx;

    //            if(m_passingVertices[i]) {
    //                std::cout << "[Passing: " << i << ", " << distances[i] <<
    //                "]" << std::endl;
    //            } else {
    //                std::cout << "[Not Passed: " << i << ", " << distances[i]
    //                << "]" << std::endl;
    //            }
  }

  // update promising vertices
  m_promisingVertices.reserve(num_verts);
  m_promisingVertices.resize(0);

  for (int i = 0; i < num_verts; ++i) {
    // mark as promising for next step capturing
    if (distances[i] < -0.5 * dx) m_promisingVertices.push_back(i);
  }
}

void StrandDynamicTraits::updateInterfaceSegments() {
  auto controller = m_strand.getParent()->getFluidScriptingController();
  if (!controller) return;

  const int num_verts = m_strand.getNumVertices();
  const Scalar dx = controller->getDx();

  for (int i = 0; i < num_verts; ++i) {
    const Vec3x pos = m_strand.getVertex(i);
    const Scalar phi = controller->get_dense_liquid_signed_distance(pos);

    m_outsideVertices[i] = phi >= -0.5 * dx;

    // we use different criteria for vertices near to bulk liquid
    m_nearVertices[i] = phi <= 0.0;
  }

  const int num_edges = m_strand.getNumEdges();
  m_insideEdges.reserve(num_edges);
  m_insideEdges.resize(0);

  for (int i = 0; i < num_edges; ++i) {
    if (!(m_outsideVertices[i] && m_outsideVertices[i + 1])) {
      m_insideEdges.push_back(i);
    }

    m_interfacingEdges[i] = m_outsideVertices[i] ^ m_outsideVertices[i + 1];
  }
}

}  // namespace strandsim
