/**
 * \copyright 2014 Danny Kaufman
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSIM_STRANDDYNAMICTRAITS_HH
#define STRANDSIM_STRANDDYNAMICTRAITS_HH

#include "../Core/BandMatrixFwd.hh"
#include "../Core/Definitions.hh"
#include "../Core/ElasticStrandParameters.hh"

namespace strandsim {

class ElasticStrand;
class DOFScriptingController;
class StrandRenderer;

class StrandDynamicTraits {
 public:
  StrandDynamicTraits(ElasticStrand& strand);

  ~StrandDynamicTraits();

  void resizeSelf();

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {}

  const std::vector<int>& getInsideEdges() const { return m_insideEdges; }

  bool isEdgeInterfacing(int eidx) const { return m_interfacingEdges[eidx]; }

  bool isVertexNear(int vidx) const { return m_nearVertices[vidx]; }

  bool isVertexInside(int vidx) const { return !m_outsideVertices[vidx]; }

  bool isVertexPassing(int vidx) const { return m_passingVertices[vidx]; }

  const std::vector<int>& getPromisingVertices() const {
    return m_promisingVertices;
  }

  const VecXx& getDisplacements() const { return m_displacements; }

  VecXx& getDisplacements() { return m_displacements; }
  Vec3x getDisplacement(IndexType vtx) const {
    return m_displacements.segment<3>(4 * vtx);
  }
  void setDisplacement(IndexType vtx, const Vec3x& disp) {
    m_displacements.segment<3>(4 * vtx) = disp;
  }
  VecXx& getInitialDisplacements() { return m_initialDisplacements; }
  VecXx& getAccelerations() { return m_accelerations; }
  const VecXx& getAccelerations() const { return m_accelerations; }
  Vec3x getAcceleration(IndexType vtx) const {
    return m_accelerations.segment<3>(4 * vtx);
  }

  void invalidatePhysics() {
    m_futureForcesUpToDate = false;
    m_futureJacobianUpToDate = false;
    m_DOFmassesUpToDate = false;
    m_flowMassesUpToDate = false;
  }
  void invalidateFuturePhysics() {
    m_futureForcesUpToDate = false;
    m_futureJacobianUpToDate = false;
  }

  // Debug drawing
  int flashForRendering();
  void addHitPoint(const Vec3x& point);
  void clearDebugDrawing();

  // Flags

  bool isFutureJacobianUpToDate() const { return m_futureJacobianUpToDate; }

  void setFutureJacobianUpToDate(bool ok) {
    this->m_futureJacobianUpToDate = ok;
  }

  bool isFutureForcesUpToDate() const { return m_futureForcesUpToDate; }

  void updateInterfaceSegments();
  void checkPassingVertices();

  // Dynamic
  void computeDOFMasses();
  void computeFlowMasses();
  void computeViscousForceCoefficients(Scalar dt);
  void computeFutureJacobian(bool withStretch = false, bool withViscous = true,
                             bool butOnlyForBendingModes = false,
                             bool dump_data = false,
                             std::ostream& dump_stream = std::cout);
  void computeLHS(Scalar dt, bool withStretch, bool withViscous,
                  bool dump_data = false,
                  std::ostream& dump_stream = std::cout);
  void computeFutureForces(bool withStretch = false, bool withViscous = true,
                           bool butOnlyForBendingModes = false,
                           bool dump_data = false,
                           std::ostream& dump_stream = std::cout);
  void computeFutureConservativeEnergy(bool withStretch = false);
  const VecXx& getDOFMasses() const;
  const VecXx& getFlowMasses() const;
  void addMassMatrixTo(JacobianMatrixType& J, Scalar multiplier = 1.0) const;
  void multiplyByMassMatrix(VecXx& F, Scalar multiplier = 1.0) const;
  void multiplyByFlowMatrix(VecXx& F) const;

  void acceptGuess();
  void acceptDisplacements();

  void setScriptingController(DOFScriptingController* controller) {
    m_scriptingController = controller;
  }

  DOFScriptingController* getScriptingController() {
    return m_scriptingController;
  }

  bool isImmune(int edge) const;
  void setImmune(int edge, bool immune) { m_immuneEdges[edge] = immune; }

  /*
   * @brief Checks if current state has NaNs; it that case replace the step with
   * rigid motion.
   *
   * This is supposed to be called after a step has left the pre-step position
   * in m_futureState.
   */
  void nanFailSafe();

 private:
  ElasticStrand& m_strand;

  // Displacements only needed for current state so we put it here. Another
  // solution is to make as shared pointer as for the Jacobian
  VecXx m_displacements;
  VecXx m_accelerations;

  // Only used in DynamicStepper/NLCImpulse
  VecXx m_initialDisplacements;

  // This could actually replace fixed vertices
  DOFScriptingController* m_scriptingController;

  VecXx m_DOFmasses;
  VecXx m_flowMasses;

  // Flags
  bool m_futureJacobianUpToDate;
  bool m_futureForcesUpToDate;
  bool m_DOFmassesUpToDate;
  bool m_flowMassesUpToDate;

  std::vector<bool> m_immuneEdges;
  std::vector<bool> m_outsideVertices;
  std::vector<bool> m_nearVertices;
  std::vector<int> m_insideEdges;
  std::vector<int> m_promisingVertices;
  std::vector<bool> m_interfacingEdges;
  std::vector<bool> m_passingVertices;

 public:
  StrandRenderer* m_renderer;  // DEBUG

  static bool s_freeMemoryAfterComputeJacobian;
};

}  // namespace strandsim

#endif  // STRANDSIM_STRANDDYNAMICTRAITS_HH
