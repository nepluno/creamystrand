/**
 * \copyright 2011 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "ElasticStrand.hh"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

#include "../Collision/CollisionUtils.hh"
#include "../Collision/OrientedBoundingBox.hh"
#include "../Dependencies/ReferenceFrames.hh"
#include "../Dynamic/DOFScriptingController.hh"
#include "../Dynamic/FluidScriptingController.hh"
#include "../Dynamic/ImplicitStepper.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Dynamic/StrandImplicitManager.hh"
#include "../Forces/BendingForce.hh"
#include "../Forces/ForceAccumulator.hh"
#include "../Forces/GravitationForce.hh"
#include "../Forces/StretchingForce.hh"
#include "../Forces/TwistingForce.hh"
#include "../Utils/EigenSerialization.hh"
#include "../Utils/MathUtilities.hh"
#include "../Utils/TextLog.hh"
#include "ElasticStrandUtils.hh"
#include "LinearSolver.hh"

namespace strandsim {

ElasticStrand::ElasticStrand(const VecXx& dofs, ParametersType& parameters,
                             CollisionParameters& collision_parameters,
                             DOFScriptingController* controller,
                             int globalIndex, const Vec3x& initRefFrame1)
    : m_globalIndex(globalIndex),           //
      m_numVertices((int)dofs.size() / 4),  //
      m_parameters(parameters),             //
      m_collisionParameters(collision_parameters),
      m_currentState(new StrandState(dofs, m_parameters.getBendingMatrixBase(),
                                     m_parameters)),  //
      m_futureState(new StrandState(dofs, m_parameters.getBendingMatrixBase(),
                                    m_parameters)),  //
      m_savedState(new StrandState(dofs, m_parameters.getBendingMatrixBase(),
                                   m_parameters)),
      m_dynamics(NULL),
      m_parent(NULL),
      m_stepper(NULL),
      m_requiresExactJacobian(false),  //
      m_projectJacobian(false),
      m_activelySimulated(true),      //
      m_isClumpCenterLine(false)      //
{
  // Allocate the Jacobian matrix and store it in a shared pointer.
  m_currentState->m_totalJacobian = m_futureState->m_totalJacobian =
      std::shared_ptr<JacobianMatrixType>(new JacobianMatrixType);

  if (controller) {
    createDynamics();
    m_dynamics->setScriptingController(controller);
  } else {
    exit(-1);
  }

  resizeInternals();

  freezeRestShape(
      0, m_numEdges);  // for now the rest shape is the shape in which the
                       // strand is created, unless modified later on.

  m_flow_reservoir.setZero();
}

ElasticStrand::ElasticStrand(const VecXx& dofs, const VecXx& area_dofs,
                             ParametersType& parameters,
                             CollisionParameters& collision_parameters,
                             DOFScriptingController* controller,
                             int globalIndex, const Vec3x& initRefFrame1)
    : m_globalIndex(globalIndex),           //
      m_numVertices((int)dofs.size() / 4),  //
      m_parameters(parameters),             //
      m_collisionParameters(collision_parameters),
      m_currentState(new StrandState(dofs, area_dofs,
                                     m_parameters.getBendingMatrixBase(),
                                     m_parameters)),  //
      m_savedState(new StrandState(dofs, area_dofs,
                                   m_parameters.getBendingMatrixBase(),
                                   m_parameters)),  //
      m_futureState(new StrandState(dofs, area_dofs,
                                    m_parameters.getBendingMatrixBase(),
                                    m_parameters)),  //
      m_dynamics(NULL),
      m_parent(NULL),
      m_stepper(NULL),
      m_requiresExactJacobian(false),  //
      m_projectJacobian(false),
      m_activelySimulated(true),      //
      m_isClumpCenterLine(false)      //
{
  // Allocate the Jacobian matrix and store it in a shared pointer.
  m_currentState->m_totalJacobian = m_futureState->m_totalJacobian =
      std::shared_ptr<JacobianMatrixType>(new JacobianMatrixType);

  if (controller) {
    createDynamics();
    m_dynamics->setScriptingController(controller);
  } else {
    exit(-1);
  }

  resizeInternals();

  freezeRestShape(
      0, m_numEdges);  // for now the rest shape is the shape in which the
                       // strand is created, unless modified later on.

  m_flow_reservoir.setZero();
}

ElasticStrand::~ElasticStrand() {
  // Statics/Dynamic should be deleted before StrandStates
  // ( as they may have stuff stored on them, such as StaticSollisionsInfo )

  delete m_dynamics;

  delete m_currentState;
  delete m_futureState;

  // Forces deleted either in StaticStrandManager or via shared pointers
}

void ElasticStrand::createDynamics() {
  if (!m_dynamics) {
    m_dynamics = new StrandDynamicTraits(*this);
  }
}

// To be called on creation
void ElasticStrand::resizeInternals() {
  m_currentState->resizeSelf();
  m_futureState->resizeSelf();

  const IndexType ndofs = getCurrentDegreesOfFreedom().size();
  m_numVertices = (ndofs + 1) / 4;
  m_numEdges =
      m_numVertices -
      1;  // This assumes open topology; allowing for loops would start here...

  if (m_dynamics) {
    m_dynamics->resizeSelf();
  }

  m_restLengths.resize(m_numEdges);
  m_restKappas.resize(m_numEdges);
  m_restTwists.resize(m_numEdges);
  m_vertexMasses.resize(m_numVertices);
  m_VoronoiLengths.resize(m_numVertices);
  m_invVoronoiLengths.resize(m_numVertices);
  m_restNeighborDistances.resize(m_numVertices);
}

void ElasticStrand::invalidatePhysics() {
  if (m_dynamics) m_dynamics->invalidatePhysics();
}

void ElasticStrand::invalidateFuturePhysics() {
  if (m_dynamics) m_dynamics->invalidateFuturePhysics();
}
void ElasticStrand::invalidateCurrentGeometry() {}

bool ElasticStrand::isThetaGoaled(int vtx) const {
  if (!m_dynamics) return false;

  return m_dynamics->getScriptingController()->isThetaGoaled(vtx);
}

bool ElasticStrand::isVertexGoaled(int vtx) const {
  if (!m_dynamics) return false;

  return m_dynamics->getScriptingController()->isVertexGoaled(vtx);
}

bool ElasticStrand::isVertexFreezed(int vtx) const {
  if (!m_dynamics) return false;

  return m_dynamics->getScriptingController()->isVertexFreezed(vtx);
}

bool ElasticStrand::isThetaFreezed(int vtx) const {
  if (!m_dynamics) return false;

  return m_dynamics->getScriptingController()->isThetaFreezed(vtx);
}

// Take the current geometry as rest shape
void ElasticStrand::freezeRestShape(unsigned begin, unsigned end,
                                    Scalar damping) {
  // Fix rest lengths
  for (IndexType vtx = begin; vtx < end; ++vtx)
    m_restLengths[vtx] = (1. - damping) * m_currentState->m_lengths[vtx] +
                         damping * m_restLengths[vtx];
  updateEverythingThatDependsOnRestLengths();

  for (IndexType vtx = begin; vtx < end; ++vtx) {
    m_restKappas[vtx] = (1. - damping) * m_currentState->m_kappas[vtx] +
                        damping * m_restKappas[vtx];
    m_restTwists[vtx] = (1. - damping) * m_currentState->m_twists[vtx] +
                        damping * m_restTwists[vtx];
  }

  invalidatePhysics();
}

// Set rest shape from dofs
void ElasticStrand::setRestShape(const VecXx& dofs, unsigned begin,
                                 unsigned end, Scalar damping) {
  VecXx backup = getCurrentDegreesOfFreedom();
  m_currentState->setDegreesOfFreedom(dofs);

  const Vec3x& initRefFrame = m_currentState->getReferenceFrame1(0);
  m_currentState->m_referenceFrames1.storeInitialFrames(initRefFrame);

  freezeRestShape(begin, end, damping);

  for (IndexType vtx = begin + 1; vtx < end; ++vtx) {
    m_restTwists[vtx] = m_restTwists[vtx - 1] +
                        clamp2Pi(m_restTwists[vtx] - m_restTwists[vtx - 1]);
  }

  m_currentState->setDegreesOfFreedom(backup);

  invalidateCurrentGeometry();
  invalidatePhysics();
}

Scalar ElasticStrand::getContactAngle() const {
  const Scalar default_val = 40.8 / 180.0 * M_PI;
  if (!m_parent) return default_val;

  auto controller = m_parent->getFluidScriptingController();
  if (!controller) return default_val;

  return controller->getContactAngle();
}

Scalar ElasticStrand::getFlowYieldStress(const VecXx& color) const {
  const Scalar default_val = 0.0;
  if (!m_parent) return default_val;

  auto controller = m_parent->getFluidScriptingController();
  if (!controller) return default_val;

  return controller->getYieldStress(color);
}

Scalar ElasticStrand::getFlowBehaviorIndex(const VecXx& color) const {
  const Scalar default_val = 1.0;
  if (!m_parent) return default_val;

  auto controller = m_parent->getFluidScriptingController();
  if (!controller) return default_val;

  return controller->getFlowBehaviorIndex(color);
}

Scalar ElasticStrand::getFlowConsistencyIndex(const VecXx& color) const {
  const Scalar default_val = 8.9e-3;
  if (!m_parent) return default_val;

  auto controller = m_parent->getFluidScriptingController();
  if (!controller) return default_val;

  return controller->getViscosity(color);
}

void ElasticStrand::updateEverythingThatDependsOnRestLengths() {
  // Total rest length
  m_totalRestLength = 0.0;
  for (IndexType vtx = 0; vtx < m_numEdges; ++vtx)
    m_totalRestLength += m_restLengths[vtx];

  for (IndexType vtx = 0; vtx < m_numVertices; ++vtx) {
    if (vtx == 0 || vtx == m_numVertices - 1)
      m_restNeighborDistances[vtx] = 0.;
    else {
      const Scalar a = m_restLengths[vtx - 1];
      const Scalar b = m_restLengths[vtx];
      const Scalar cost = cos(m_parameters.getMinBendingAngle());
      m_restNeighborDistances[vtx] = sqrt(a * a + b * b - 2.0 * a * b * cost);
    }
  }

  // Compute Voronoi lengths
  m_VoronoiLengths[0] = 0.5 * m_restLengths[0];
  for (IndexType vtx = 1; vtx < m_numEdges; ++vtx)
    m_VoronoiLengths[vtx] = 0.5 * (m_restLengths[vtx - 1] + m_restLengths[vtx]);
  m_VoronoiLengths[m_numEdges] = 0.5 * m_restLengths[m_numVertices - 2];

  // Compute masses and inverse of Voronoi lengths
  for (IndexType vtx = 0; vtx < m_numVertices; ++vtx) {
    m_vertexMasses[vtx] = m_parameters.getDensity() * m_VoronoiLengths[vtx] *
                          M_PI * m_parameters.getRadiusA(vtx, m_numVertices) *
                          m_parameters.getRadiusB(vtx, m_numVertices);

    m_invVoronoiLengths[vtx] = 1.0 / m_VoronoiLengths[vtx];
  }

  invalidatePhysics();
}

Scalar ElasticStrand::getFutureSurfaceFlowMassAtEdge(IndexType i) const {
  const auto controller = m_parent->getFluidScriptingController();
  if (!controller) return 0.0;

  const int num_components = controller->getNumComponents();
  VecXx color = m_stepper->flowNewComponents().segment(i * num_components,
                                                       num_components) +
                m_stepper->flowNewComponents().segment((i + 1) * num_components,
                                                       num_components);

  make_gibbs_simplex(color);

  const Scalar density = controller->getDensity(color);
  return density * getFutureSurfaceFlowVolumeAtEdge(i);
}

Scalar ElasticStrand::getCurrentSurfaceFlowMassAtEdge(IndexType i) const {
  const auto controller = m_parent->getFluidScriptingController();
  if (!controller) return 0.0;

  const int num_components = controller->getNumComponents();
  VecXx color = m_stepper->flowNewComponents().segment(i * num_components,
                                                       num_components) +
                m_stepper->flowNewComponents().segment((i + 1) * num_components,
                                                       num_components);

  make_gibbs_simplex(color);

  const Scalar density = controller->getDensity(color);
  return density * getCurrentSurfaceFlowVolumeAtEdge(i);
}

Scalar ElasticStrand::getFutureSurfaceFlowVolumeAtEdge(IndexType i) const {
  if (!m_parent) return 0.0;

  const Scalar area = (m_futureState->m_area_dofs.getVertex(i) +
                       m_futureState->m_area_dofs.getVertex(i + 1)) *
                      0.5;
  return area * m_restLengths[i];
}

Scalar ElasticStrand::getCurrentSurfaceFlowVolumeAtEdge(IndexType i) const {
  if (!m_parent) return 0.0;

  const Scalar area = (m_currentState->m_area_dofs.getVertex(i) +
                       m_currentState->m_area_dofs.getVertex(i + 1)) *
                      0.5;
  return area * m_restLengths[i];
}

Scalar ElasticStrand::getCurrentSurfaceFlowMass(const IndexType i) const {
  if (!m_parent) return 0.0;

  const auto controller = m_parent->getFluidScriptingController();
  if (!controller) return 0.0;

  const int num_components = controller->getNumComponents();
  const VecXx& color =
      m_stepper->flowComponents().segment(i * num_components, num_components);

  const Scalar density = controller->getDensity(color);
  const Scalar area = m_currentState->m_area_dofs.getVertex(i);
  const Scalar mass = density * area * m_VoronoiLengths[i];

  return mass;
}

Scalar ElasticStrand::getFutureSurfaceFlowMass(const IndexType i) const {
  if (!m_parent) return 0.0;

  const auto controller = m_parent->getFluidScriptingController();
  if (!controller) return 0.0;

  const int num_components = controller->getNumComponents();
  const VecXx& color = m_stepper->flowNewComponents().segment(
      i * num_components, num_components);

  const Scalar density = controller->getDensity(color);
  const Scalar area = m_futureState->m_area_dofs.getVertex(i);
  const Scalar mass = density * area * m_VoronoiLengths[i];

  return mass;
}

void ElasticStrand::setEdgeRestLength(const IndexType vtx,
                                      const Scalar newrestlength) {
  assert(vtx < m_numEdges);

  m_restLengths[vtx] = newrestlength;

  if (0 == vtx) {
    // If we change the rest length of a fixed edge, move the vertices
    setVertex(1, getVertex(0) + newrestlength * getEdgeVector(0).normalized());
  }
}

VecXx ElasticStrand::getFlowNewComponents(int pidx) const {
  const int num_components =
      m_stepper->flowNewComponents().size() / m_numVertices;

  return m_stepper->flowNewComponents().segment(pidx * num_components,
                                                num_components);
}

void ElasticStrand::setEdgesRestLength(const Scalar newRestLength) {
  for (IndexType vtx = 0; vtx < m_numEdges; ++vtx)
    setEdgeRestLength(vtx, newRestLength);
  updateEverythingThatDependsOnRestLengths();

  invalidatePhysics();
}

void ElasticStrand::setRadius(const Scalar radius_a, const Scalar radius_b) {
  m_parameters.setRadii(radius_a, radius_b);

  for (IndexType vtx = 0; vtx < m_numVertices; ++vtx) {
    m_vertexMasses[vtx] = m_parameters.getDensity() * m_VoronoiLengths[vtx] *
                          M_PI * m_parameters.getRadiusA(vtx, m_numVertices) *
                          m_parameters.getRadiusB(vtx, m_numVertices);
  }

  invalidatePhysics();
}

void ElasticStrand::setStiffness(const Scalar youngs) {
  m_parameters.setYoungsModulus(youngs);

  invalidatePhysics();
}

void ElasticStrand::ackParametersChanged() {
  updateEverythingThatDependsOnRestLengths();
}

void ElasticStrand::setParameters(const ElasticStrandParameters& parameters) {
  m_parameters = parameters;

  ackParametersChanged();
}

void ElasticStrand::setParameters(double i_radiusA, double i_radiusB,
                                  double i_rootRM, double i_tipRM,
                                  double i_youngsModulus, double i_shearModulus,
                                  double i_density, double i_viscosity,
                                  double i_airDrag) {
  m_parameters.setRadii(i_radiusA, i_radiusB);
  m_parameters.setYoungsModulus(i_youngsModulus);
  m_parameters.setShearModulus(i_shearModulus);
  m_parameters.setViscosity(i_viscosity);
  m_parameters.setDensity(i_density);
  m_parameters.setAirDrag(i_airDrag);
  m_parameters.setRadiusMultiplier(0, i_rootRM);
  m_parameters.setRadiusMultiplier(1, i_tipRM);

  ackParametersChanged();
}

void ElasticStrand::setFutureDegreesOfFreedom(const VecXx& dof) {
  getFutureState().setDegreesOfFreedom(dof);

  invalidateFuturePhysics();
}

void ElasticStrand::setFutureAreaDegreesOfFreedom(const VecXx& dof) {
  getFutureState().setAreaDegreesOfFreedom(dof);

  invalidateFuturePhysics();
}

void ElasticStrand::setCurrentAreaDegreesOfFreedom(const VecXx& dof) {
  getCurrentState().setAreaDegreesOfFreedom(dof);

  invalidatePhysics();
}

void ElasticStrand::setSavedDegreesOfFreedom(const VecXx& dof) {
  getSavedState().setDegreesOfFreedom(dof);
}

void ElasticStrand::setCurrentDegreesOfFreedom(const VecXx& dof) {
  getCurrentState().setDegreesOfFreedom(dof);

  invalidateCurrentGeometry();
}

Scalar ElasticStrand::getUnsignedAngleToMajorRadius(int vtx,
                                                    const Vec3x& vec) const {
  if (vtx + 1 == m_numVertices) --vtx;

  const Vec3x& edge = getEdgeVector(vtx).normalized();
  const Vec3x& orth = (vec - vec.dot(edge) * edge);
  const Scalar north = orth.norm();
  if (isSmall(north)) return 0.;
  return std::acos(clamp(getMaterialFrame2(vtx).dot(orth / north), -1., 1.));
}

Scalar ElasticStrand::getSignedAngleToMajorRadius(int vtx,
                                                  const Vec3x& vec) const {
  if (vtx + 1 == m_numVertices) --vtx;

  const Vec3x& edge = getEdgeVector(vtx).normalized();
  const Vec3x& orth = (vec - vec.dot(edge) * edge);

  const Scalar cosa = getMaterialFrame2(vtx).dot(orth);
  const Scalar sina = -getMaterialFrame1(vtx).dot(orth);
  return std::atan2(sina, cosa);
}

// Add ForceT's force on theta only to the VecXx
template <typename ForceT>
void ElasticStrand::accumulateEFThetaOnly(Scalar& thetaE, VecXx& thetaF,
                                          StrandState& geometry) const {
  typename ForceT::LocalThetaForceType localF;

  for (IndexType vtx = ForceT::s_first; vtx < m_numVertices - ForceT::s_last;
       ++vtx) {
    thetaE += ForceT::localEnergy(*this, geometry, vtx);

    ForceT::computeLocal(localF, *this, geometry, vtx);
    ForceT::addInPosition(thetaF, vtx, localF);
  }
}

// Add ForceT's Jacobian on theta only to the tridiagonal matrix
template <typename ForceT>
void ElasticStrand::accumulateJThetaOnly(TriDiagonalMatrixType& thetaJ,
                                         StrandState& geometry) const {
  typename ForceT::LocalThetaJacobianType localJ;

  for (IndexType vtx = ForceT::s_first; vtx < m_numVertices - ForceT::s_last;
       ++vtx) {
    ForceT::computeLocal(localJ, *this, geometry, vtx);
    ForceT::addInPosition(thetaJ, vtx, localJ);
  }
}

void ElasticStrand::accumulateRuntimeForces(RuntimeForceBase::Quantities q,
                                            StrandState& geometry) {
  for (auto force = m_runtimeForces.begin(); force != m_runtimeForces.end();
       ++force) {
    if (*force) {
      (*force)->accumulate(q, geometry, *this);
    }
  }

  for (auto force = m_sharedRuntimeForces.begin();
       force != m_sharedRuntimeForces.end(); ++force) {
    if (*force) {
      (*force)->accumulate(q, geometry, *this);
    }
  }
}

Scalar ElasticStrand::getCurrentTotalLength() const {
  Scalar totalLength = 0.0;
  for (IndexType vtx = 0; vtx < m_numEdges; ++vtx) {
    totalLength += m_currentState->m_lengths[vtx];
  }

  return totalLength;
}

Scalar ElasticStrand::getFutureTotalLength() const {
  Scalar totalLength = 0.0;
  for (IndexType vtx = 0; vtx < m_numVertices - 1; ++vtx) {
    totalLength += m_futureState->m_lengths[vtx];
  }

  return totalLength;
}

void ElasticStrand::findBendElasticLimitMax(Scalar& bendElasticLimit) const {
  // find the max BendElastic limit per strand
  for (IndexType vtx = 1; vtx < m_numEdges; ++vtx) {
    for (int i = 0; i < 4; ++i) {
      bendElasticLimit = std::max(bendElasticLimit,
                                  fabs(Scalar(m_currentState->m_kappas[vtx][i] -
                                              m_restKappas[vtx][i])));
    }
  }
}

void ElasticStrand::applyPlasticDeformation(Scalar stretchElasticLimit,
                                            Scalar bendElasticLimit,
                                            Scalar twistElasticLimit) {
  // Rest lengths plastic deformation and everything that depend on them: total
  // rest length, Voronoi lengths, masses...
  for (int vtx = 0; vtx < m_numEdges; ++vtx)
    if (m_currentState->m_lengths[vtx] >
        m_restLengths[vtx] + stretchElasticLimit)
      m_restLengths[vtx] = m_currentState->m_lengths[vtx] - stretchElasticLimit;
  updateEverythingThatDependsOnRestLengths();

  for (IndexType vtx = 0; vtx < m_numEdges; ++vtx) {
    for (int i = 0; i < 4; ++i) {
      if (m_currentState->m_kappas[vtx][i] >
          m_restKappas[vtx][i] + bendElasticLimit) {
        m_restKappas[vtx][i] =
            m_currentState->m_kappas[vtx][i] - bendElasticLimit;
      }
      if (m_currentState->m_kappas[vtx][i] <
          m_restKappas[vtx][i] - bendElasticLimit) {
        m_restKappas[vtx][i] =
            m_currentState->m_kappas[vtx][i] + bendElasticLimit;
      }
    }

    if (m_currentState->m_twists[vtx] > m_restTwists[vtx] + twistElasticLimit) {
      m_restTwists[vtx] = m_currentState->m_twists[vtx] - twistElasticLimit;
    }
    if (m_currentState->m_twists[vtx] < m_restTwists[vtx] - twistElasticLimit) {
      m_restTwists[vtx] = m_currentState->m_twists[vtx] + twistElasticLimit;
    }
  }

  invalidatePhysics();
}

Vec3x ElasticStrand::getGoalVertex(IndexType vtx) const {
  if (!m_dynamics) return Vec3x::Zero();

  return m_dynamics->getScriptingController()->getVertexGoal(vtx);
}

Scalar ElasticStrand::getGoalTheta(IndexType vtx) const {
  if (!m_dynamics) return 0.0;

  return m_dynamics->getScriptingController()->getThetaGoal(vtx);
}

void ElasticStrand::filterFutureGeometryByRestLength(const double epsilon,
                                                     bool allowCompression) {
  const Scalar ca = std::cos(m_parameters.getMinBendingAngle());

  Vec3x xaN;
  Vec3x xaP;

  for (int i = 0; i < m_numVertices; ++i) {
    // ignore freezed vertices
    if (isVertexFreezed(i) || isVertexGoaled(i)) {
      xaN = xaP = m_futureState->getVertex(i);
      continue;
    }

    Vec3x xbN = m_futureState->getVertex(i);
    Vec3x xbNrev = xbN;

    const Scalar lP = getEdgeRestLength(i - 1);

    if (i >= 2) {
      // first check if the angle is correct
      const Vec3x& xaP2 = m_futureState->getVertex(i - 2);
      const Scalar lP2 = getEdgeRestLength(i - 2);

      const Scalar lJ = (xbNrev - xaP2).norm();
      const Scalar lJmin = sqrt(lP * lP + lP2 * lP2 - 2. * lP * lP2 * ca);
      if (lJ > 0. && lJ < lJmin) {
        Scalar ratio2 = lJmin / lJ;
        xbNrev = xaP2 + (xbNrev - xaP2) * ratio2;
      }
    }

    const Scalar lN = (xbNrev - xaP).norm();
    Scalar ratio = 1.;

    if (!isSmall(lN)) {
      if (lN > lP + epsilon) {
        ratio = (lP + epsilon) / lN;
      } else if (!allowCompression && lN < lP - epsilon) {
        ratio = (lP - epsilon) / lN;
      }
    }

    // compute and store revised delta
    xbNrev = xaN + (xbNrev - xaP) * ratio;

    m_futureState->setVertex(i, xbNrev);

    xaP = xbN;
    xaN = xbNrev;
  }

  invalidateFuturePhysics();
}

void ElasticStrand::addRuntimeForce(RuntimeForceBase* force) {
  m_runtimeForces.push_back(force);

  invalidatePhysics();
}

Scalar ElasticStrand::getDt() const { return m_parent->getDt(); }

void ElasticStrand::addSharedRuntimeForce(
    const std::shared_ptr<RuntimeForceBase> force) {
  m_sharedRuntimeForces.push_back(force);

  invalidatePhysics();
}

void ElasticStrand::clearSharedRuntimeForces() {
  m_sharedRuntimeForces.clear();
}

Scalar ElasticStrand::getCurvilinearAbscissa(int vtx,
                                             Scalar localAbscissa) const {
  Scalar s = 0;

  for (auto i = 0; i < vtx; ++i) {
    s += getEdgeRestLength(i);
  }
  if (vtx + 1 < m_numVertices) {
    s += localAbscissa * getEdgeRestLength(vtx);
  }

  return s;
}

void ElasticStrand::getLocalAbscissa(const Scalar curvilinearAbscissa, int& vtx,
                                     Scalar& localAbscissa) const {
  for (localAbscissa = curvilinearAbscissa, vtx = 0;
       vtx + 1 < m_numVertices && m_restLengths[vtx] < localAbscissa;
       localAbscissa -= m_restLengths[vtx++])
    ;

  if (vtx + 1 == m_numVertices) {
    localAbscissa = 0;
  } else {
    localAbscissa /= m_restLengths[vtx];
  }
}

void ElasticStrand::getAABB(StrandState& geometry, unsigned elementID,
                            Vec3x& min, Vec3x& max,
                            CollisionParameters::CollisionType type) {
  static const Vec3x unit = Vec3x::Ones();

  const std::pair<Vec3x, Vec3x>& aaBB = geometry.getAABB(elementID);
  const Scalar rad =
      m_collisionParameters.collisionsRadius(type, elementID, m_numVertices);
  const Scalar flow_area =
      std::max(geometry.getAreaDegreesOfFreedom(elementID),
               geometry.getAreaDegreesOfFreedom(elementID + 1));
  const Scalar max_dist =
      (1. + 0.5 * m_collisionParameters.m_cohesionTheta) * sqrt(flow_area);

  min = aaBB.first - (rad + max_dist) * unit;
  max = aaBB.second + (rad + max_dist) * unit;
}

void ElasticStrand::getAABB(unsigned elementID, Vec3x& min, Vec3x& max,
                            CollisionParameters::CollisionType type) {
  getAABB(getCurrentState(), elementID, min, max, type);
}

void ElasticStrand::getFutureAABB(unsigned elementID, Vec3x& min, Vec3x& max,
                                  CollisionParameters::CollisionType type) {
  getAABB(getFutureState(), elementID, min, max, type);
}

void ElasticStrand::getSegment(unsigned elementID, Vec3x& start,
                               Vec3x& end) const {
  assert(elementID + 1 < getNumVertices());
  start = getVertex(elementID);
  end = getVertex(elementID + 1);
}

template <class RuntimeForceT>
void ElasticStrand::removeSharedForcesOfType() {
  unsigned nDeleted = 0;

  for (std::list<std::shared_ptr<RuntimeForceBase> >::iterator force =
           m_sharedRuntimeForces.begin();
       force != m_sharedRuntimeForces.end();) {
    if (std::dynamic_pointer_cast<RuntimeForceT>(*force)
            .use_count())  // If the casted shared_ptr is not empty
    {
      // No need to delete a shared_ptr!
      force = m_sharedRuntimeForces.erase(force);
      ++nDeleted;
    } else
      force++;
  }

  if (nDeleted) invalidatePhysics();
}

void ElasticStrand::addClumpingAttractor(const ElasticStrand* strand) {
  m_clumpingAttractors.push_back(strand);

  invalidatePhysics();
}

void ElasticStrand::removeClumpingAttractor(const ElasticStrand* strand) {
  ClumpingAttractorsIteratorType whereInAttractorContainer = std::find(
      m_clumpingAttractors.begin(), m_clumpingAttractors.end(), strand);
  m_clumpingAttractors.erase(whereInAttractorContainer);

  invalidatePhysics();
}

void ElasticStrand::setClumpingAttractors(
    ClumpingAttractorsContainerType& attractors) {
  m_clumpingAttractors = attractors;

  invalidatePhysics();
}

void ElasticStrand::removeAllClumpingAttractors() {
  m_clumpingAttractors.clear();

  invalidatePhysics();
}

void ElasticStrand::removeRuntimeForce(const RuntimeForceBase* toDelete) {
  for (auto force = m_runtimeForces.begin(); force != m_runtimeForces.end();) {
    if (toDelete == *force) {
      force = m_runtimeForces.erase(force);
      // Already deleted by StaticStrandManager
    } else
      force++;
  }

  invalidatePhysics();
}

void ElasticStrand::removeSharedRuntimeForce(
    const std::shared_ptr<RuntimeForceBase> toDelete) {
  for (RuntimeForcesIteratorType force =
           find(m_sharedRuntimeForces.begin(), m_sharedRuntimeForces.end(),
                toDelete);
       force != m_sharedRuntimeForces.end();
       force = find(force, m_sharedRuntimeForces.end(), toDelete)) {
    force = m_sharedRuntimeForces.erase(force);
  }
}

template <class RuntimeForceT>
void ElasticStrand::removeForcesOfType() {
  for (std::list<RuntimeForceBase*>::iterator force = m_runtimeForces.begin();
       force != m_runtimeForces.end();)
    if (dynamic_cast<const RuntimeForceT*>(*force)) {
      delete *force;
      force = m_runtimeForces.erase(force);
    } else
      force++;

  invalidatePhysics();
}

void ElasticStrand::printCurrentEnergies() {
  const VecXx oldF = m_currentState->m_totalForce;
  const Scalar oldE = m_currentState->m_totalEnergy;

  std::cout << " Total: " << m_currentState->m_totalEnergy
            << " ; |F| = " << m_currentState->m_totalForce.norm() << std::endl;
  std::cout << " --- > " << oldF.norm() / getNumVertices() << " / "
            << oldF.lpNorm<Eigen::Infinity>() << std::endl;

  m_currentState->m_totalEnergy = 0.0;
  m_currentState->m_totalForce.setZero();
  accumulateEF<StretchingForce<> >(*m_currentState);
  std::cout << " Stretching: " << m_currentState->m_totalEnergy
            << " ; |F| = " << m_currentState->m_totalForce.norm() << std::endl;

  m_currentState->m_totalEnergy = 0.0;
  m_currentState->m_totalForce.setZero();
  accumulateEF<TwistingForce<> >(*m_currentState);
  std::cout << " Twisting: " << m_currentState->m_totalEnergy
            << " ; |F| = " << m_currentState->m_totalForce.norm() << std::endl;

  m_currentState->m_totalEnergy = 0.0;
  m_currentState->m_totalForce.setZero();
  accumulateEF<BendingForce<> >(*m_currentState);
  std::cout << " Bending: " << m_currentState->m_totalEnergy
            << " ; |F| = " << m_currentState->m_totalForce.norm() << std::endl;

  m_currentState->m_totalEnergy = 0.0;
  m_currentState->m_totalForce.setZero();
  accumulateEF<GravitationForce>(*m_currentState);
  std::cout << " Grav: " << m_currentState->m_totalEnergy
            << " ; |F| = " << m_currentState->m_totalForce.norm() << std::endl;

  std::cout << " Kappas " << std::endl;
  for (unsigned i = 0; i < m_numEdges; ++i) {
    std::cout << m_currentState->m_kappas.get()[i];
  }
  std::cout << std::endl << " RestKappas " << std::endl;
  for (unsigned i = 0; i < m_numEdges; ++i) {
    std::cout << m_restKappas[i];
  }
  std::cout << std::endl;
  std::cout << " Twists " << std::endl;
  for (unsigned i = 0; i < m_numEdges; ++i) {
    std::cout << " " << m_currentState->m_twists.get()[i];
  }
  std::cout << std::endl << " RestTwists " << std::endl;
  for (unsigned i = 0; i < m_numEdges; ++i) {
    std::cout << " " << m_restTwists[i];
  }
  std::cout << std::endl;

  m_currentState->m_totalForce = oldF;
  m_currentState->m_totalEnergy = oldE;
}

std::ostream& operator<<(std::ostream& os, const ElasticStrand& strand) {
  const VecXx& dofs = strand.getCurrentDegreesOfFreedom();
  os << '{';
  for (int i = 0; i < strand.m_numEdges; i++) {
    os << '{' << dofs[4 * i] << ", " << dofs[4 * i + 1] << ", "
       << dofs[4 * i + 2] << "}, ";
    // os << strand.m_currentState->m_degreesOfFreedom[4 * i + 3] << ', ';
  }
  os << '{' << dofs[4 * (strand.m_numEdges)] << ", "
     << dofs[4 * (strand.m_numEdges) + 1] << ", "
     << dofs[4 * (strand.m_numEdges) + 2] << '}';
  os << '}';

  return os;
}

/**
 * This will be called as a safeguard if relaxThetas fails: put all the
 * non-fixed thetas to zero. If someone has a better idea...
 */
void ElasticStrand::resetThetas() {
  m_futureState->setThetas(VecXx(m_numEdges));
}

bool ElasticStrand::serializeTo(std::ostream& os) const {
  std::cerr << "ElasticStrand::serializeTo not supported temporarily"
            << std::endl;

  return false;
}

bool ElasticStrand::deserializeFrom(std::istream& is) {
  std::cerr << "ElasticStrand::deserializeFrom not supported temporarily"
            << std::endl;

  return false;
}

}  // namespace strandsim
