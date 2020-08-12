/**
 * \copyright 2014 Danny Kaufman, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSIM_STRANDIMPLICITMANAGER_HH
#define STRANDSIM_STRANDIMPLICITMANAGER_HH

#include <list>
#include <memory>
#include <set>
#include <vector>

#include "../../bogus/Interfaces/MecheEigenInterface.hpp"
#include "../Collision/ProximityCollision.hh"
#include "../Core/Definitions.hh"
#include "../Utils/SpatialHashMapFwd.hh"
#include "CohesionTable.hh"
#include "SimulationParameters.hh"

namespace strandsim {

class ElasticStrand;
class ImplicitStepper;
class DOFScriptingController;
class CollisionDetector;
class ElementProxy;
class MeshScriptingController;
class FluidScriptingController;
class ConstraintScriptingController;
class CollisionBase;
class LevelSetForce;

class StrandImplicitManager {
 public:
  //! Simulation step timings -- in milliseconds
  struct SubStepCallback {
    virtual void executeCallback() = 0;
    virtual void projectConstraint(const std::vector<ImplicitStepper*>&) = 0;
  };

  struct SubStepTimings {
    SubStepTimings()
        : prepare(0),
          hairHairCollisions(0),
          dynamics(0),
          meshHairCollisions(0),
          processCollisions(0),
          solve(0) {}

    double prepare;             // Executing controllers, etc
    double hairHairCollisions;  // Hair hair collision detection
    double dynamics;            // Assembling + pre-solving linear systems
    double meshHairCollisions;  // Mesh hair collision detection
    double processCollisions;   // Updating collision gradients
    double solve;               // Proper solve

    double sum() const {
      return prepare + hairHairCollisions + dynamics + meshHairCollisions +
             processCollisions + solve;
    }
  };

  //! collision detection timings -- in milliseconds
  struct CDTimings {
    CDTimings()
        : buildBVH(0.),
          findCollisionsBVH(0.),
          narrowPhase(0.),
          upateHashMap(0.),
          processHashMap(0.) {}

    double buildBVH;
    double findCollisionsBVH;
    double narrowPhase;
    double upateHashMap;
    double processHashMap;

    void reset() {
      buildBVH = 0.;
      findCollisionsBVH = 0.;
      narrowPhase = 0.;
      upateHashMap = 0.;
      processHashMap = 0.;
    }
  };

  //! Statistics for on simulation timestep
  struct SolverStats {
    unsigned totConstraints;
    unsigned maxConstraints;
    unsigned maxObjects;
    double maxError;
    double maxTime;

    void reset();
  };

  typedef std::vector<SubStepTimings> StepTimings;
  typedef std::list<StepTimings> Timings;

  StrandImplicitManager(
      const std::vector<ElasticStrand*>& strands,
      const std::map<std::pair<int, int>, std::set<std::pair<int, int> > >&
          collision_free,
      const std::vector<std::shared_ptr<MeshScriptingController> >&
          meshScripting_controllers,
      const std::vector<std::shared_ptr<FluidScriptingController> >&
          fluidScripting_controllers,
      const std::vector<ConstraintScriptingController*>&
          constraintScripting_controllers,
      Scalar startTime, Scalar dt, const SimulationParameters& params,
      SubStepCallback* sub_callback);

  virtual ~StrandImplicitManager();

  //! Step the simulation forward by m_dt
  /*! \param substepsInDt the number of simulation substeps to perform */
  void execute(int total_num_substeps, int total_substep_id,
               const Scalar total_substep_dt);

  //! Performs a simulation substep
  void step(int total_num_substeps, int total_substep_id);

  Scalar getCFL() const;

  void setDt(Scalar dt) { m_dt = dt; }
  Scalar getDt() const { return m_dt; }
  Scalar getTime() const { return m_time; }
  const std::string& getOutputDirectory() const { return m_output_directory; }

  void setOutputDirectory(const std::string& dir) { m_output_directory = dir; }

  void setTime(Scalar m_time) { this->m_time = m_time; }
  void updateParameters(const SimulationParameters& params);

  void print(const StepTimings& timings) const;
  void print(const Timings& timings) const;

  const Timings& timings() const { return m_timings; }

  void drawContacts() const;

  const std::vector<ElasticStrand*>& getStrands() const { return m_strands; }

  std::shared_ptr<FluidScriptingController> getFluidScriptingController(
      int idx = 0) const {
    if (idx >= m_fluidScriptingControllers.size()) {
      return NULL;
    } else {
      return m_fluidScriptingControllers[idx];
    }
  }

 private:
  template <typename StreamT>
  void print(const SubStepTimings& timings) const;
  template <typename StreamT>
  void print(const SolverStats& stats) const;

  void printMemStats();

  //! Map between a index in the simulation to an index in a colliding groyup
  typedef std::map<unsigned, unsigned> IndicesMap;
  //! Colliding group: set of strands and contacts that should be solved
  //! together
  typedef std::pair<IndicesMap, ProximityCollisions> CollidingGroup;

  //! Prepares the substep, executes the externals objects controllers
  void step_prepare(Scalar dt);
  //! Solves the unconstrained dynamics of each strand
  void step_dynamics(int total_num_substeps, int total_substep_id, Scalar dt);
  //! Rewind and redo the unconstrained dynamics of each strand (possibly with
  //! new forces)
  void redo_step_dynamics(int total_num_substeps, int total_substep_id,
                          Scalar dt);
  //! Analyses the current collisions and compute the colliding groups
  void step_processCollisions(Scalar dt);

  //! Solves each flow friction group
  void step_solveFlowFrictions(int total_num_substeps, int total_substep_id);

  //! Solves each colliding group
  void step_solveCollisions(int total_num_substeps, int total_substep_id);

  //! Returns wether a strand needs to be solved using bogus.
  /*! Will be true if the strand is subject to at least one contact or hard
   * constraint */
  bool needsExternalSolve(unsigned strandIdx) const;

  //! Returns wether a strand needs to be solved using bogus.
  /*! Will be true if the strand is subject to at least one contact or hard
   * constraint */
  bool needsElasticExternalSolve(unsigned strandIdx) const;

  //! Setup the contacts and constraints for a colliding group
  bool assembleBogusFrictionProblem(
      CollidingGroup& collisionGroup, bogus::MecheFrictionProblem& mecheProblem,
      std::vector<ProximityCollisions>& externalContacts,
      std::vector<unsigned>& globalIds,
      std::vector<ProximityCollision*>& colPointers, VecXx& vels,
      VecXx& worldImpulses, VecXx& impulses, VecXx& adhesions, VecXx& filters,
      VecXu& startDofs, VecXu& nDofs, int& numSubSys,
      bool herschelBulkleyProblem, bool ignoreMutualCollision);
  //! Cleanup a friction problem and updates the strands with the new velocities
  //! if \p accept is true
  int postProcessBogusFrictionProblem(
      bool updateVelocity, CollidingGroup& collisionGroup,
      const bogus::MecheFrictionProblem& mecheProblem,
      const std::vector<unsigned>& globalIds,
      const std::vector<ProximityCollision*>& colPointers, VecXx& vels,
      VecXx& worldImpulses, VecXx& impulses, VecXu& startDofs, VecXu& nDofs,
      std::vector<Scalar>& newton_residuals, int total_num_substeps,
      int total_substep_id);
  //! Proper solving of the MecheFrictionProblem
  Scalar solveBogusFrictionProblem(bogus::MecheFrictionProblem& mecheProblem,
                                   const std::vector<unsigned>& globalIds,
                                   bool asFailSafe, bool herschelBulkleyProblem,
                                   bool doFrictionShrinking, VecXx& vels,
                                   VecXx& worldImpulses, VecXx& impulses,
                                   int& numSubSys);

  //! Solve the contacts and constraints on a single object
  Scalar solveSingleObject(std::vector<ProximityCollisions>& externalContacts,
                           const unsigned objectIdx, bool asFailSafe,
                           bool herschelBulkleyProblem, bool updateVelocity,
                           bool ignoreMutualCollision,
                           std::vector<Scalar>& newtonResiduals,
                           int& numNewtonIters, int total_num_substeps,
                           int total_substep_id);
  //! Solve the contacts and constraints on a colliding group
  Scalar solveCollidingGroup(CollidingGroup& cg,
                             std::vector<ProximityCollisions>& externalContacts,
                             bool asFailSafe, bool herschelBulkleyProblem,
                             bool updateVelocity, bool ignoreMutualCollision,
                             std::vector<Scalar>& newtonResiduals,
                             int& numNewtonIters, int total_num_substeps,
                             int total_substep_id);

  void residualStats(const std::string& name,
                     const std::vector<Scalar>& residuals);

  //! Mesh/hair collision detection
  void setupMeshHairCollisions(Scalar dt);
  //! Hair/hair collision detection
  void setupHairHairCollisions(Scalar dt);

  //! Proximity mesh/hair collision detection
  void doProximityMeshHairDetection(Scalar dt);
  //! Continuous-time mesh/hair collisision detection
  void doContinuousTimeDetection(Scalar dt);

  bool isCollisionInvariantCT(Scalar dt);

  //! Adds an external contact on strand \p strIdx, edge \p edgeIdx, abscissa \p
  //! abscissa
  /*! \return whether this collision has been accepted */
  bool addExternalContact(const unsigned strIdx, const unsigned edgeIdx,
                          const Scalar abscissa,
                          const ProximityCollision& collision);

  //! Computes the deformation gradient of a strand at one contact point, ie
  //! dq/dx
  void computeDeformationGradient(ProximityCollision::Object& object) const;
  //! Setup the local frame for one contact and calls
  //! computeDeformationGradient() for each object
  void setupDeformationBasis(ProximityCollision& collision) const;

  //! Transform a mutual collision into an external contact on the
  //! (onFirstObject ? first : second) object
  void makeExternalContact(ProximityCollision& c, bool onFirstObject,
                           std::vector<ProximityCollisions>& contacts);
  void makeExternalContact(ProximityCollision& c, bool onFirstObject);
  void makeElasticExternalContact(ProximityCollision& c, bool onFirstObject);
  //! Discards mesh/hair collisions that are unlikely to be activated
  void pruneExternalCollisions(std::vector<ProximityCollisions>& contacts);

  //! Discards hair/hair collisions that are unlikely to be activated
  void pruneCollisions(const ProximityCollisions& origMutualCollisions,
                       ProximityCollisions& mutualCollisions,
                       const Scalar stochasticPruning, bool elastic);
  //! Computes the colliding groups using a graph walking algorithm
  void computeCollidingGroups(const ProximityCollisions& mutualCollisions,
                              Scalar dt, bool elastic);

  void exportStrandRestShapes(const std::string& fileName) const;

  Scalar m_time;  //!< Current time
  Scalar m_dt;    //!< Time per "frame"
  SimulationParameters m_params;
  std::string m_output_directory;

  const std::vector<ElasticStrand*>& m_strands;
  const std::map<std::pair<int, int>, std::set<std::pair<int, int> > >&
      m_collision_free;
  std::vector<ImplicitStepper*> m_steppers;
  const std::vector<std::shared_ptr<MeshScriptingController> >&
      m_meshScriptingControllers;
  const std::vector<std::shared_ptr<FluidScriptingController> >&
      m_fluidScriptingControllers;
  std::vector<std::shared_ptr<LevelSetForce> > m_levelSetForces;

  const std::vector<ConstraintScriptingController*>&
      m_constraintScriptingControllers;
  //    std::vector<TriangularMesh*> m_triangularMeshes;
  std::vector<ElementProxy*>
      m_elementProxies;  //!< List of all proxies that should be inserted in the
                         //!< BVH
  CollisionDetector* m_collisionDetector;  //!< BVH-based collision detector

  std::vector<ProximityCollisions>
      m_externalContacts;  //!< External contacts on each strand
  std::vector<ProximityCollisions>
      m_elasticExternalContacts;  //!< External contacts under elastic rule only

  ProximityCollisions m_mutualContacts;         //!< List of all mutual contacts
  ProximityCollisions m_elasticMutualContacts;  //!< List of all mutual contacts
                                                //!< under elastic rule only
  //! Structure for storing collisions and forces, useful for drawaing and
  //! warm-starting solver
  ProximityCollisionDatabase m_collisionDatabase;

  std::vector<CollidingGroup> m_collidingGroups;
  std::vector<CollidingGroup> m_elasticCollidingGroups;

  std::vector<unsigned> m_globalIds;

  //! Index of colliding group in which each strand should be. Can be -1.
  std::vector<int> m_collidingGroupsIdx;
  std::vector<int> m_elasticCollidingGroupsIdx;
  //! Set of strands that share a common bilat constraint
  std::vector<std::set<unsigned> > m_commonBilateralConstraints;

  //!< Spatial Hash Map for hair/hair proximity collision detetection
  typedef SpatialHashMap<ElasticStrand, unsigned, true> SpatialHashMapT;
  SpatialHashMapT* m_hashMap;

  SubStepCallback* m_substep_callback;
  // Stat gathering
  Timings m_timings;
  SolverStats m_stats;

  Scalar m_mem_usage_accu;
  Scalar m_mem_usage_divisor;
  void printProblemStats(const StrandImplicitManager::SubStepTimings& timings);
  unsigned m_statExternalContacts, m_statMutualCollisions, m_statTotalSubsteps;
  unsigned m_num_nonlinear_iters, m_num_contact_solves, m_max_nonlinear_iters,
      m_max_perstep_nonlinear_iters;
  int m_num_ct_hair_hair_col;
  int m_unconstrained_NewtonItrs;
  CDTimings m_cdTimings;
};

StrandImplicitManager::SubStepTimings operator+(
    const StrandImplicitManager::SubStepTimings& lhs,
    const StrandImplicitManager::SubStepTimings& rhs);
StrandImplicitManager::SubStepTimings operator/(
    const StrandImplicitManager::SubStepTimings& lhs, const Scalar rhs);

}  // namespace strandsim

#endif  // STRANDSIM_STRANDIMPLICITMANAGER_HH
