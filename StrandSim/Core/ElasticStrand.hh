/**
 * \copyright 2011 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ELASTICSTRAND_HH_
#define ELASTICSTRAND_HH_

#include <list>
#include <memory>
#include <set>

#include "../Forces/ForceBase.hh"
#include "../Utils/ThreadUtils.hh"
#include "BandMatrixFwd.hh"
#include "CollisionParameters.hh"
#include "ElasticStrandParameters.hh"
#include "StrandBase.hh"
#include "StrandState.hh"

#ifndef M_PI
#define M_PI 3.141592653589793238462650288
#endif

namespace strandsim {

template <typename ViscousT>
class StretchingForce;
template <typename ViscousT>
class BendingForce;
template <typename ViscousT>
class TwistingForce;
template <typename ViscousT>
class FixingForce;
template <typename ForceT>
class ForceAccumulator;

class GravitationForce;
class AirDragForce;
class MassDampingForce;
class StrandDynamicTraits;
class DOFScriptingController;
class StrandImplicitManager;
class ImplicitStepper;

struct ProximityCollision;

class ElasticStrand : public StrandBase {
 public:
  typedef ElasticStrandParameters ParametersType;
  typedef std::list<std::shared_ptr<RuntimeForceBase> >
      RuntimeForcesContainerType;
  typedef RuntimeForcesContainerType::iterator RuntimeForcesIteratorType;
  typedef std::list<const ElasticStrand*> ClumpingAttractorsContainerType;
  typedef std::list<const ElasticStrand*>::iterator
      ClumpingAttractorsIteratorType;
  typedef std::list<const ElasticStrand*>::const_iterator
      ClumpingAttractorsConstIteratorType;

  ElasticStrand(const VecXx& dofs, ParametersType& parameters,
                CollisionParameters& collision_parameters,
                DOFScriptingController* controller, int globalIndex = -1,
                const Vec3x& initRefFrame1 = Vec3x());

  ElasticStrand(const VecXx& dofs, const VecXx& area_dofs,
                ParametersType& parameters,
                CollisionParameters& collision_parameters,
                DOFScriptingController* controller, int globalIndex = -1,
                const Vec3x& initRefFrame1 = Vec3x());

  virtual ~ElasticStrand();

  IndexType getNumVertices() const { return m_numVertices; }
  IndexType getNumEdges() const { return m_numEdges; }

  const std::vector<ProximityCollision>& getPenaltyContacts() const;
  // Acknowledges that current physics are not valid anymore. Which happens more
  // often than one may think. This could be due to several things, for instance
  // a change in physical parameters or external forces. For statics, it means
  // that we may no longer be in an equilibrium positions ; we have to resume
  // updating
  void invalidatePhysics();
  // Invalidate physics of future state only
  void invalidateFuturePhysics();
  // Acknowledge that the strand's shape has changed
  void invalidateCurrentGeometry();

  VecXx getFlowNewComponents(int pidx) const;

  ////////////////////////////////////////////////////////////////////////////////
  // Access to current state
  const VecXx& getCurrentAreaDegreesOfFreedom() const {
    return m_currentState->m_area_dofs.get();
  }
  const VecXx& getCurrentDegreesOfFreedom() const {
    return m_currentState->m_dofs.get();
  }
  const VecXx& getSavedDegreesOfFreedom() const {
    return m_savedState->m_dofs.get();
  }
  void setCurrentDegreesOfFreedom(const VecXx& dof);
  void setSavedDegreesOfFreedom(const VecXx& dof);

  template <typename Derived>
  void getCurrentVertices(Eigen::MatrixBase<Derived>& vertices) const {
    typedef typename Eigen::internal::traits<Derived>::Scalar OutScalar;

    const VecXx& dofs = getCurrentDegreesOfFreedom();
    for (IndexType vtx = 0; vtx < m_numVertices; ++vtx) {
      vertices.template segment<3>(3 * vtx) =
          dofs.segment<3>(4 * vtx).cast<OutScalar>();
    }
  }
  Vec2x& getReservoir() { return m_flow_reservoir; }
  const Vec2x& getReservoir() const { return m_flow_reservoir; }
  Vec3x getGoalVertex(IndexType vtx) const;
  Scalar getGoalTheta(IndexType vtx) const;

  Vec3x getVertex(IndexType vtx) const {
    return m_currentState->getVertex(vtx);
  }
  void setVertex(IndexType vtx, const Vec3x& point) {
    m_currentState->setVertex(vtx, point);
  }
  Scalar getTheta(const IndexType vtx) const {
    assert(vtx < m_numEdges);
    return m_currentState->getTheta(vtx);
  }

  // Flow Behaviors
  Scalar getContactAngle() const;
  Scalar getFlowYieldStress(const VecXx& color) const;
  Scalar getFlowBehaviorIndex(const VecXx& color) const;
  Scalar getFlowConsistencyIndex(const VecXx& color) const;

  void setTheta(const IndexType vtx, const Scalar newTheta) {
    m_currentState->setTheta(vtx, newTheta);
    invalidatePhysics();
  }

  Scalar getVoronoiLength(const IndexType vtx) const {
    return m_VoronoiLengths[vtx];
  }

  Scalar getCurrentDynamicVoronoiLength(const IndexType vtx) const {
    return m_currentState->m_dynamic_voronoi_lengths[vtx];
  }

  Scalar getFutureDynamicVoronoiLength(const IndexType vtx) const {
    return m_futureState->m_dynamic_voronoi_lengths[vtx];
  }

  const Vec3x& getCurrentTangent(int vtx) const {
    return m_currentState->m_tangents[vtx];
  }
  const Vec3xArray& getCurrentReferenceFrames1() {
    return m_currentState->m_referenceFrames1.get();
  }
  const Vec3xArray& getCurrentReferenceFrames2() {
    return m_currentState->m_referenceFrames2.get();
  }
  const Vec3xArray& getCurrentMaterialFrames1() {
    return m_currentState->m_materialFrames1.get();
  }
  const Vec3xArray& getCurrentMaterialFrames2() {
    return m_currentState->m_materialFrames2.get();
  }
  const std::vector<Scalar>& getCurrentReferenceTwists() const {
    return m_currentState->m_referenceTwists.get();
  }
  void setCurrentReferenceFrames1(const Vec3xArray& reff1) {
    m_currentState->m_referenceFrames1.set(reff1);
  }
  void setCurrentReferenceFrames2(const Vec3xArray& reff2) {
    m_currentState->m_referenceFrames2.set(reff2);
  }
  void setCurrentReferenceTwists(const std::vector<Scalar>& reft) {
    m_currentState->m_referenceTwists.set(reft);
  }

  Vec3x getMaterialFrame1(int vtx) const {
    return m_currentState->getMaterialFrame1(vtx);
  }

  Vec3x getMaterialFrame2(int vtx) const {
    return m_currentState->getMaterialFrame2(vtx);
  }

  Vec3x getEdgeCenter(int edge) const {
    int next_vert = std::min(edge + 1, getNumVertices() - 1);
    return (getVertex(edge) + getVertex(next_vert)) * 0.5;
  }

  Vec3x getEdgeDirectionAtVertex(int vtx) const {
    const int n = getNumVertices();
    if (vtx == 0)
      return getEdgeVector(0).normalized();
    else if (vtx == n - 1)
      return getEdgeVector(vtx - 1).normalized();
    else
      return (getEdgeVector(vtx) + getEdgeVector(vtx - 1)).normalized();
  }

  Vec3x getEdgeVector(int vtx) const {
    return m_currentState->getEdgeVector(vtx);
  }

  Scalar getFutureFlowDOFArea(int vtx) const {
    return m_futureState->getAreaDegreesOfFreedom()(vtx);
  }

  Scalar getCurrentFlowDOFArea(int vtx) const {
    return m_currentState->getAreaDegreesOfFreedom()(vtx);
  }

  void setCurrentFlowDOFArea(int vtx, const Scalar& area) {
    m_currentState->setAreaDegreesOfFreedom(vtx, area);
  }

  Scalar getFutureFlowHeight(int vtx) const {
    return m_futureState->m_height[vtx];
  }

  Scalar getCurrentFlowHeight(int vtx) const {
    return m_currentState->m_height[vtx];
  }

  Scalar getFutureFlowLaplaceHeight(int vtx) const {
    return m_futureState->m_laplace_height[vtx];
  }

  Vec3x getFutureEdgeVector(int vtx) const {
    return m_futureState->getEdgeVector(vtx);
  }

  Vec3x getFutureEdgeDirection(int vtx) const {
    return m_futureState->getEdgeVector(vtx) /
           m_futureState->getEdgeLength(vtx);
  }

  bool isVertexFreezed(int vtx) const;
  bool isThetaFreezed(int vtx) const;
  bool isVertexGoaled(int vtx) const;
  bool isThetaGoaled(int vtx) const;

  Scalar getUnsignedAngleToMajorRadius(int vtx, const Vec3x& vec) const;
  Scalar getSignedAngleToMajorRadius(int vtx, const Vec3x& vec) const;

  Scalar getCurrentTotalLength() const;
  void freezeRestShape(unsigned begin, unsigned end, Scalar damping = 0.);

  ///////////////////////////////////////////
  // Future geometry

  Vec3x getFutureVertex(IndexType vtx) const {
    return m_futureState->getVertex(vtx);
  }

  const VecXx& getFutureDegreesOfFreedom() const {
    return m_futureState->m_dofs.get();
  }
  void setFutureDegreesOfFreedom(const VecXx& dof);
  void setFutureAreaDegreesOfFreedom(const VecXx& dof);
  void setCurrentAreaDegreesOfFreedom(const VecXx& dof);

  Scalar getFutureTotalLength() const;

  void filterFutureGeometryByRestLength(const double epsilon = 0.,
                                        bool allowCompression = false);

  Scalar getTotalEnergy() const { return m_currentState->m_totalEnergy; }

  const VecXx& getTotalForces() const { return m_currentState->m_totalForce; }

  VecXx& getTotalForces() { return m_currentState->m_totalForce; }

  VecXx& getFutureTotalForces() { return m_futureState->m_totalForce; }

  const JacobianMatrixType& getTotalJacobian() const {
    return *(m_currentState->m_totalJacobian);
  }

  JacobianMatrixType& getTotalJacobian() {
    return *(m_currentState->m_totalJacobian);
  }

  Scalar getNewTotalEnergy() const { return m_futureState->m_totalEnergy; }

  const VecXx& getNewTotalForces() const { return m_futureState->m_totalForce; }

  VecXx& getNewTotalForces() { return m_futureState->m_totalForce; }

  void swapStates() {
    // std::cout << "SWAP [" << m_globalIndex << "]" <<  std::endl;
    std::swap(m_currentState, m_futureState);
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Changing the shape

  void findBendElasticLimitMax(Scalar& bendElasticLimit) const;

  void applyPlasticDeformation(Scalar stretchElasticLimit,
                               Scalar bendElasticLimit,
                               Scalar twistElasticLimit);

  void updateEverythingThatDependsOnRestLengths();

  ////////////////////////////////////////////////////////////////////////////////
  // Changing the forces

  void addRuntimeForce(RuntimeForceBase* force);
  void addSharedRuntimeForce(const std::shared_ptr<RuntimeForceBase> force);
  void clearSharedRuntimeForces();

  const std::list<RuntimeForceBase*>& getRuntimeForces() const {
    return m_runtimeForces;
  }

  const RuntimeForcesContainerType& getSharedRuntimeForces() const {
    return m_sharedRuntimeForces;
  }

  void addClumpingAttractor(const ElasticStrand* strand);
  void removeClumpingAttractor(const ElasticStrand* strand);
  void setClumpingAttractors(ClumpingAttractorsContainerType& attractors);
  void removeAllClumpingAttractors();
  size_t numberOfClumpingAttractors() { return m_clumpingAttractors.size(); }

  void removeRuntimeForce(const RuntimeForceBase* toDelete);
  void removeSharedRuntimeForce(
      const std::shared_ptr<RuntimeForceBase> toDelete);

  template <class RuntimeForceT>
  void removeForcesOfType();

  template <class RuntimeForceT>
  void removeSharedForcesOfType();

  void setIsClumpCenterLine(bool isit) { m_isClumpCenterLine = isit; }

  bool getIsClumpCenterLine() const { return m_isClumpCenterLine; }

  /////////////////////////////
  // Rest Shape and physical parameters

  Scalar getEdgeRestLength(const IndexType vtx) const {
    return m_restLengths[std::max((IndexType)0,
                                  std::min(vtx, (IndexType)(m_numEdges - 1)))];
  }

  void setEdgeRestLength(const IndexType vtx, const Scalar newrestlength);

  void setEdgesRestLength(const Scalar newRestLength);

  void setRestShape(const VecXx& dofs, unsigned begin, unsigned end,
                    Scalar damping = 0.);

  Scalar getRadiusA(int vtx) const {
    return m_parameters.getRadiusA(vtx, getNumVertices());
  }

  Scalar getRadiusB(int vtx) const {
    return m_parameters.getRadiusB(vtx, getNumVertices());
  }

  void setRadius(const Scalar radius_a, const Scalar radius_b);

  void multiplyRadii(Scalar m) {
    m_parameters.internalRadiusMultiplier() *= m;

    invalidatePhysics();
    invalidateCurrentGeometry();
  }

  Scalar getRestTwist(const IndexType vtx) const { return m_restTwists[vtx]; }

  void setRestTwist(const IndexType vtx, const Scalar restTwist) {
    m_restTwists[vtx] = restTwist;
    invalidatePhysics();
  }

  Vec4x getKappaBar(const IndexType vtx) const { return m_restKappas[vtx]; }

  void setKappaBar(const IndexType vtx, const Vec4x& kappaBar) {
    m_restKappas[vtx] = kappaBar;
    invalidatePhysics();
  }

  Scalar getStiffness() const { return m_parameters.getYoungsModulus(); }

  void setStiffness(const Scalar youngs);

  void setParameters(double i_radiusA, double i_radiusB, double i_rootRM,
                     double i_tipRM, double i_youngsModulus,
                     double i_shearModulus, double i_density,
                     double i_viscosity, double i_airDrag);

  void setParameters(const ElasticStrandParameters& parameters);

  const ElasticStrandParameters& getParameters() const { return m_parameters; }
  ElasticStrandParameters& getParameters() { return m_parameters; }

  //! Informs the strand that its parameters have been modified
  void ackParametersChanged();

  const CollisionParameters& collisionParameters() const {
    return m_collisionParameters;
  }

  Scalar getTotalRestLength() const { return m_totalRestLength; }

  Scalar getVertexMass(IndexType i) const { return m_vertexMasses[i]; }

  Scalar getDt() const;

  Scalar getCurrentSurfaceFlowMass(IndexType i) const;

  Scalar getFutureSurfaceFlowMass(IndexType i) const;

  Scalar getFutureSurfaceFlowMassAtEdge(IndexType i) const;

  Scalar getCurrentSurfaceFlowMassAtEdge(IndexType i) const;

  Scalar getFutureSurfaceFlowVolumeAtEdge(IndexType i) const;

  Scalar getCurrentSurfaceFlowVolumeAtEdge(IndexType i) const;

  Scalar getEdgeFlowInertia(IndexType i) const {
    const Scalar a = m_parameters.getRadiusA(i, getNumVertices());
    const Scalar b = m_parameters.getRadiusB(i, getNumVertices());
    const Scalar mass1 =
        (getFutureSurfaceFlowMass(i) + getFutureSurfaceFlowMass(i + 1)) * 0.5;

    const Scalar ha00 = a + getFutureFlowHeight(i);
    const Scalar ha01 = b + getFutureFlowHeight(i);
    const Scalar ha10 = a + getFutureFlowHeight(i + 1);
    const Scalar ha11 = b + getFutureFlowHeight(i + 1);

    return 0.25 * mass1 *
           (square(a) + square(b) +
            (square(ha00) + square(ha01) + square(ha10) + square(ha11)) * 0.5);
  }

  Scalar getCurrentEdgeFlowInertia(IndexType i) const {
    const Scalar a = m_parameters.getRadiusA(i, getNumVertices());
    const Scalar b = m_parameters.getRadiusB(i, getNumVertices());
    const Scalar mass1 =
        (getCurrentSurfaceFlowMass(i) + getCurrentSurfaceFlowMass(i + 1)) * 0.5;

    const Scalar ha00 = a + getCurrentFlowHeight(i);
    const Scalar ha01 = b + getCurrentFlowHeight(i);
    const Scalar ha10 = a + getCurrentFlowHeight(i + 1);
    const Scalar ha11 = b + getCurrentFlowHeight(i + 1);

    return 0.25 * mass1 *
           (square(a) + square(b) +
            (square(ha00) + square(ha01) + square(ha10) + square(ha11)) * 0.5);
  }

  Scalar getEdgeInertia(IndexType i) const {
    const Scalar a = m_parameters.getRadiusA(i, getNumVertices());
    const Scalar b = m_parameters.getRadiusB(i, getNumVertices());
    const Scalar mass0 =
        m_parameters.getDensity() * M_PI * a * b * m_restLengths[i];
    const Scalar mass1 =
        (getFutureSurfaceFlowMass(i) + getFutureSurfaceFlowMass(i + 1)) * 0.5;

    const Scalar ha00 = a + getFutureFlowHeight(i);
    const Scalar ha01 = b + getFutureFlowHeight(i);
    const Scalar ha10 = a + getFutureFlowHeight(i + 1);
    const Scalar ha11 = b + getFutureFlowHeight(i + 1);

    return 0.25 * (mass0 * (square(a) + square(b)) +
                   mass1 * (square(a) + square(b) +
                            (square(ha00) + square(ha01) + square(ha10) +
                             square(ha11)) *
                                0.5));
  }

  Scalar getCurvilinearAbscissa(int vtx, Scalar localAbscissa) const;
  void getLocalAbscissa(const Scalar curvilinearAbscissa, int& vtx,
                        Scalar& localAbscissa) const;

  //////////////////////////////////////////////
  // States, static, dynamics

  StrandState& getSavedState() { return *m_savedState; }

  const StrandState& getSavedState() const { return *m_savedState; }

  StrandState& getCurrentState() { return *m_currentState; }

  const StrandState& getCurrentState() const { return *m_currentState; }

  const StrandState& getFutureState() const { return *m_futureState; }

  StrandState& getFutureState() { return *m_futureState; }

  bool canDoDynamics() const { return m_dynamics; }
  void createDynamics();
  const StrandDynamicTraits& dynamics() const { return *m_dynamics; }
  StrandDynamicTraits& dynamics() { return *m_dynamics; }

  bool requiresExactJacobian() const { return m_requiresExactJacobian; }
  void requireExactJacobian(bool b) { m_requiresExactJacobian = b; }

  bool projectsJacobian() const { return m_projectJacobian; }
  void projectJacobian(bool b) { m_projectJacobian = b; }

  bool activelySimulated() const { return m_activelySimulated; }

  void setActivelySimulated(bool simulated) { m_activelySimulated = simulated; }

  void printCurrentEnergies();

  // Operates on future state's thetas
  void resetThetas();

  ///////////////
  // Segments and bounding boxes

  void getAABB(
      unsigned elementID, Vec3x& min, Vec3x& max,
      CollisionParameters::CollisionType type = CollisionParameters::SELF);
  void getFutureAABB(
      unsigned elementID, Vec3x& min, Vec3x& max,
      CollisionParameters::CollisionType type = CollisionParameters::SELF);

  void getSegment(unsigned elementID, Vec3x& start, Vec3x& end) const;
  // Begin and end of edge iterators, used in SpatialHashMap
  unsigned subsamples_begin() const { return 0; }
  unsigned subsamples_end() const { return m_numEdges; }

  //////////////////////
  // Serialization

  bool serializeTo(std::ostream& os) const;
  bool deserializeFrom(std::istream& is);

  template <class Archive>
  void serialize(Archive& ar, const int version) {
    ar & m_currentState->m_dofs;
    ar & m_currentState->m_referenceFrames1;
    ar & m_currentState->m_referenceFrames2;
    ar & m_currentState->m_referenceTwists;
    ar & m_currentState->m_referenceFrames1.getPreviousTangents();
    ar& m_parameters.dt();
    ar& m_parameters.internalRadiusMultiplier();
  }

  MutexType& mutex() { return *m_mutex; }

  int getGlobalIndex() const { return m_globalIndex; }

  void setGlobalIndex(int globalIndex) { m_globalIndex = globalIndex; }

  const Vec4xArray& getRestKappas() const { return m_restKappas; }

  const std::vector<Scalar>& getRestLengths() const { return m_restLengths; }

  const std::vector<Scalar>& getRestTwists() const { return m_restTwists; }

  void setParent(StrandImplicitManager* parent) { m_parent = parent; }

  StrandImplicitManager* getParent() const { return m_parent; }

  void setStepper(ImplicitStepper* stepper) { m_stepper = stepper; }

  ImplicitStepper* getStepper() const { return m_stepper; }

 private:
  void getAABB(StrandState& geometry, unsigned elementID, Vec3x& min,
               Vec3x& max, CollisionParameters::CollisionType type);

  const std::list<const ElasticStrand*>& getClumpingAttractors() const {
    return m_clumpingAttractors;
  }

  void resizeInternals();

  /////////////////////////////////////
  // Convenience force accumulators

  // Add ForceT's energy to the geometry
  template <typename ForceT>
  void accumulateE(StrandState& geometry) const {
    ForceAccumulator<ForceT>::accumulate(geometry.m_totalEnergy, *this,
                                         geometry);
  }

  // Add ForceT to the geometry
  template <typename ForceT>
  void accumulateF(StrandState& geometry) const {
    ForceAccumulator<ForceT>::accumulate(geometry.m_totalForce, *this,
                                         geometry);
  }

  // Add ForceT's Jacobian to the geometry
  template <typename ForceT>
  void accumulateJ(StrandState& geometry) const {
    ForceAccumulator<ForceT>::accumulate(*geometry.m_totalJacobian, *this,
                                         geometry);
  }

  // Add ForceT's energy and force to the geometry
  template <typename ForceT>
  void accumulateEF(StrandState& geometry) const {
    accumulateE<ForceT>(geometry);
    accumulateF<ForceT>(geometry);
  }

  // Add ForceT's energy, force and Jacobian to the geometry
  template <typename ForceT>
  void accumulateEFJ(StrandState& geometry) const {
    accumulateE<ForceT>(geometry);
    accumulateF<ForceT>(geometry);
    accumulateJ<ForceT>(geometry);
  }
  // Add ForceT's energy to the geometry
  template <typename ForceT>
  void accumulate(ForceBase::Quantities q, StrandState& geometry) const {
    switch (q) {
      case ForceBase::EFJ:
        accumulateJ<ForceT>(geometry);
        /* no break */
      case ForceBase::EF:
        accumulateF<ForceT>(geometry);
        /* no break */
      case ForceBase::E:
        accumulateE<ForceT>(geometry);
        break;
      case ForceBase::F:
        accumulateF<ForceT>(geometry);
        break;
      case ForceBase::J:
        accumulateJ<ForceT>(geometry);
        break;
      case ForceBase::NONE:
        break;
    };
  }

  template <typename ForceT>
  void accumulateEFThetaOnly(Scalar& thetaEnergy, VecXx& thetaForce,
                             StrandState& geometry) const;
  template <typename ForceT>
  void accumulateJThetaOnly(TriDiagonalMatrixType& J,
                            StrandState& geometry) const;

  void accumulateRuntimeForces(RuntimeForceBase::Quantities q,
                               StrandState& geometry);

  //////////////////////////////////////////////
  /**
   * Member variables
   */

  int m_globalIndex;  // Global index in the simulation

  // Size of the strand. The original belongs to m_parameters, so it can
  // correctly interpolate when asked e.g. for a variable radius.
  IndexType m_numVertices;
  IndexType m_numEdges;

  // Other physical parameters
  ParametersType& m_parameters;
  CollisionParameters& m_collisionParameters;

  // Current and future geometry. They need to be pointers to be easily swapped.
  StrandState* m_currentState;
  StrandState* m_futureState;
  StrandState* m_savedState;

  // Statics / Dynamics
  StrandDynamicTraits* m_dynamics;
  // Parent to access other strands
  StrandImplicitManager* m_parent;
  ImplicitStepper* m_stepper;

  // Rest shape
  std::vector<Scalar>
      m_restLengths;  // The following four members depend on m_restLengths,
                      // which is why updateEverythingThatDependsOnRestLengths()
                      // must be called
  std::vector<Scalar> m_restNeighborDistances;
  Scalar m_totalRestLength;
  std::vector<Scalar> m_VoronoiLengths;     // rest length around each vertex
  std::vector<Scalar> m_invVoronoiLengths;  // their inverses
  std::vector<Scalar> m_vertexMasses;
  Vec4xArray m_restKappas;
  std::vector<Scalar> m_restTwists;

  // Reservior for dripping
  Vec2x m_flow_reservoir;

  // Flags
  bool m_requiresExactJacobian;
  bool m_projectJacobian;
  bool m_activelySimulated;

  // Forces that are not built-in
  std::list<RuntimeForceBase*> m_runtimeForces;

  // Forces that are shared between strands
  RuntimeForcesContainerType m_sharedRuntimeForces;

  // Other stuff
  ClumpingAttractorsContainerType m_clumpingAttractors;

  bool m_isClumpCenterLine;
  MutexWrapper m_mutex;

  friend class Viscous;
  friend class NonViscous;
  template <typename ViscousT>
  friend class StretchingForce;
  template <typename ViscousT>
  friend class BendingForce;
  template <typename ViscousT>
  friend class TwistingForce;
  friend class GravitationForce;
  friend class AirDragForce;
  friend class MassDampingForce;
  friend class StrandDynamicTraits;
  friend class FluidDragForce;
  template <typename ViscousT>
  friend class FixingForce;

  friend std::ostream& operator<<(std::ostream& os,
                                  const ElasticStrand& strand);
};

}  // namespace strandsim

BOOST_CLASS_VERSION(strandsim::ElasticStrand, 0)

#endif /* ELASTICSTRAND_HH_ */
