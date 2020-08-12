/**
 * \copyright 2014 Danny Kaufman, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSIM_IMPLICITSTEPPER_HH
#define STRANDSIM_IMPLICITSTEPPER_HH

#include <memory>
#include <vector>

#include "../Collision/ProximityCollision.hh"
#include "../Core/StepperBase.hh"
#include "../Utils/SymmetricBandMatrixSolver.hh"

namespace strandsim {

class ElasticStrand;
struct SimulationParameters;

class ProximityCollision;
typedef std::vector<ProximityCollision> ProximityCollisions;

class ImplicitStepper {
 public:
  ImplicitStepper(ElasticStrand& strand, const SimulationParameters& params);

  virtual ~ImplicitStepper();

  //! Updates the strand degrees of freedom using the velocities stored in
  //! m_newVelocities
  /*! Checks if the strand has a high stretch energy. If this is the case and
    afterContraints if false, calls an appropriate failsafe \return whether the
    new velocities are acceptable ( stretch energy is low enough )
  */
  bool update(bool afterConstraints = false);

  void updateRodAccelerationAndVelocity();

  void updateRodAcceleration();

  void setDt(const Scalar& dt) { m_dt = dt; }

  void setFraction(const Scalar& fraction) { m_fraction = fraction; }

  //! Accept's the new strand's state
  void finalize();

  //! Scales the dynamics' linear sytem by s
  void scale(const Scalar s);

  //! Reverts the strand's state to the one at the beginning of the timestep
  void rewind();

  Scalar getCurrentFlowVelocityAtVertex(int vtx) const;

  VecXx& flowNewStrain() { return m_flow_newStrain; }

  const VecXx& flowNewStrain() const { return m_flow_newStrain; }

  VecXx& flowStrain() { return m_flow_strain; }

  const VecXx& flowStrain() const { return m_flow_strain; }

  const VecXx& getCurrentFlowVelocity() const { return m_flow_velocities; }

  const VecXx& getFutureFlowVelocity() const { return m_flow_newVelocities; }

  void setFutureFlowVelocity(int eidx, const Scalar& u_tau) {
    m_flow_newVelocities[eidx] = u_tau;
  }

  void setCurrentFlowVelocity(int eidx, const Scalar& u_tau) {
    m_flow_velocities[eidx] = u_tau;
  }

  Scalar getDt() const { return m_dt; }

  void initStepping(Scalar dt);
  //! Starts a new substep
  void startSubstep(int id, Scalar dt, Scalar fraction);
  //! Update the strand's boundary conditions controller
  void updateDofController(int subStepId, int numSubstepsPerDt);

  ElasticStrand& getStrand() { return m_strand; }

  const ElasticStrand& getStrand() const { return m_strand; }

  JacobianMatrixType& Lhs();
  const JacobianMatrixType& Lhs() const;

  const VecXx& flow_lhs() const { return m_flow_lhs; }

  const VecXx& rhs() const { return m_rhs; }

  VecXx& impulse_rhs();

  const VecXx& getDOFMasses() const;

  VecXx& flowVelocities() { return m_flow_velocities; }
  VecXx& newFlowVelocities() { return m_flow_newVelocities; }

  VecXx& velocities() { return m_velocities; }
  VecXx& newVelocities() { return m_newVelocities; }

  const VecXx& flowComponents() const { return m_flow_components; }
  const VecXx& flowNewComponents() const { return m_flow_newComponents; }

  VecXx& flowComponents() { return m_flow_components; }
  VecXx& flowNewComponents() { return m_flow_newComponents; }

  VecXx& additionalImpulses() { return m_additional_impulse; }

  const VecXx& additionalImpulses() const { return m_additional_impulse; }

  const VecXx& velocities() const { return m_velocities; }
  const VecXx& newVelocities() const { return m_newVelocities; }
  bool notSPD() const { return m_notSPD; }

  bool lastStepWasRejected() const { return m_lastStepWasRejected; }

  bool usedNonLinearSolver() { return m_usedNonlinearSolver; }

  bool refusesMutualContacts() const { return notSPD(); }

  JacobianSolver& linearSolver() { return m_linearSolver; }

  JacobianSolver& massMatrixLinearSolver() { return m_massMatrix_linearSolver; }

  //    JacobianSolver& complianceLinearSolver ()
  //    {
  //        return m_compliance_linearSolver;
  //    }

  //! Updates the current Lhs and rhs based of the m_newVelocities guess
  /*! \return whether the linear system has been updated */
  bool updateLinearSystem(const VecXx solverForces);

  unsigned numIters() { return m_newtonIter; }

  void prepareSolveNonlinear();

  bool postSolveNonlinear();

  bool prepareNewtonIteration();

  bool performNewtonIteration();
  //! Intialized future degrees of freedom, setup frames, optionally initialize
  //! length constraints
  void prepareDynamics();

  void backtrackFlowData();

  void updateAdditionalInertia();

  void recoverFromBeginning();

  void stepFlowAdvForce();

  void stepFlowData();

  void limitAdditionalImpulses(const Scalar& maxImpulses);

  void updateRHSwithImpulse(const VecXx& impulses);

  Scalar getNewtonResidual() const { return m_minErr; }

  Scalar maxAdditionalImpulseNorm(int& idx);

  Scalar getStretchMultiplier() const;
  void solveLinear();

 private:
  //! Computes linearized dynamics

  void solveLinearFlow();

  void updateFlowData();

  //! Computes the left-hand-side of the linear system of linearized dynamics at
  //! current guess
  void computeLHS(bool dump_data = false, std::ostream& = std::cout);
  //! Computes the right-hand-side of the linear system of linearized dynamics
  //! at current guess
  void computeRHS(bool dump_data = false, std::ostream& = std::cout);

  //! Computes the left-hand-side of the linear system of flow dynamics
  void computeFlowLHS();

  //! Computes the right-hand-side of the linear system of flow dynamics
  void computeFlowRHS();

  //! Copy reference frames from current state to future state
  void setupFuturesFrames();

  void clearConstraints();
  //! Allocates the length constraints
  void createConstraints();

  VecXx interpolateStrand(const Scalar& pos,
                          const std::vector<Scalar>& summed_length,
                          const VecXx& source, bool clamped,
                          const VecXx& default_val);

  Vec4x interpolateStrandVelocity(const Scalar& pos,
                                  const std::vector<Scalar>& summed_length,
                                  const VecXx& source, bool clamped,
                                  const Vec4x& default_val);

  Scalar interpolateStrand(const Scalar& pos,
                           const std::vector<Scalar>& summed_length,
                           const VecXx& source, bool clamped,
                           const Scalar& default_val);

  Scalar traceRK2Strand(const Scalar& pos, const Scalar& dt,
                        const std::vector<Scalar>& summed_length,
                        const VecXx& u_vert, bool clamped,
                        const Scalar& default_val);

  Scalar traceRK2Strand(const int ipos, const Scalar& dt,
                        const std::vector<Scalar>& summed_length,
                        const VecXx& u_vert, bool clamped,
                        const Scalar& default_val);

  void accumulateReservoir();

  void markFlowStatus();

  template <int nelem = 1, int nstep = 1>
  void mapEdgeToVertices(const VecXx& edge_vars, VecXx& vertex_vars);

  template <int nelem = 1, int nstep = 1>
  void mapVertexToEdges(const VecXx& vertex_vars, VecXx& edge_vars);

  void computeGradAtVertex(const VecXx& edge_vars, VecXx& vertex_grads);
  //! Geometric projection of the strand's future dofs to enfore edges rest
  //! lengths
  /*! \param preStep  whether this projection is done before or after dynamics
   */
  void filterGeometryLength(bool preStep);

  //! Returns the value of the stretch energy divided  by the length of the rods
  //! times its stiffness
  Scalar getLineicStretch();

  Scalar m_dt;
  Scalar m_fraction;

  const SimulationParameters& m_params;
  Scalar m_stretchingFailureThreshold;
  Scalar m_costretchResidualFailureThreshold;
  Scalar m_stretchMultiplier;

  JacobianSolver m_linearSolver;

  JacobianSolver m_massMatrix_linearSolver;  // for zeroth-order contact resolve
  JacobianMatrixType m_massMatrix;           // really just mass

  //    JacobianSolver m_compliance_linearSolver;
  //    JacobianMatrixType m_complianceMatrix;

  bool m_notSPD;
  bool m_usedNonlinearSolver;
  bool m_linearSystemIsDirty;
  bool m_lastStepWasRejected;

  VecXx m_flow_velocities;
  VecXx m_flow_newVelocities;

  VecXx m_flow_components;
  VecXx m_flow_newComponents;

  VecXx m_flow_strain;
  VecXx m_flow_newStrain;

  VecXx m_beginningVelocities;
  VecXx m_savedVelocities;
  VecXx m_velocities;
  VecXx m_newVelocities;
  VecXx m_additional_inertia;

  VecXx m_rhs;
  VecXx m_impulseRhs;  // for zeroth-order contact
  VecXx m_additional_impulse;
  VecXx m_projectionDisplacements;

  VecXx m_flow_rhs;
  VecXx m_flow_lhs;

  JacobianMatrixType m_bestLHS;
  VecXx m_bestRhs;
  VecXx m_prevRhs;

  //	std::vector<unsigned char> m_flow_status;
  //	VecXx m_flow_pressure;

  Scalar m_alpha;  // Current step length
  Scalar m_minErr;
  Scalar m_prevErr;

  Scalar m_total_flow;

 public:
  ElasticStrand& m_strand;

 private:
  unsigned m_newtonIter;
};

}  // namespace strandsim

#endif  // STRANDSIM_IMPLICITSTEPPER_HH
