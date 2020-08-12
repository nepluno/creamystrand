/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef DOFSCRIPTINGCONTROLLER_HH_
#define DOFSCRIPTINGCONTROLLER_HH_

#include <map>

#include "../Core/BandMatrixFwd.hh"

namespace strandsim {

class StrandState;

class DOFScriptingController {
 public:
  DOFScriptingController();
  DOFScriptingController(const Vec3xArray&);

  virtual ~DOFScriptingController();

  void clear() {
    m_scriptedDegreesOfFreedom.clear();
    m_scriptedGoal.clear();
  }

  bool isVertexGoaled(int vtx) const {
    return m_scriptedGoal.find(vtx * 4) != m_scriptedGoal.end();
  }

  bool isThetaGoaled(int vtx) const {
    return m_scriptedGoal.find(vtx * 4 + 3) != m_scriptedGoal.end();
  }

  bool isVertexFreezed(int vtx) const {
    return m_scriptedDegreesOfFreedom.find(vtx * 4) !=
           m_scriptedDegreesOfFreedom.end();
  }

  bool isThetaFreezed(int vtx) const {
    return m_scriptedDegreesOfFreedom.find(vtx * 4 + 3) !=
           m_scriptedDegreesOfFreedom.end();
  }

  template <int n>
  void freezeRootVertices() {
    for (int dof = 0; dof < 4 * n - 1; ++dof)
      m_scriptedDegreesOfFreedom[dof] = 0.0;
  }

  void ungoalVertices(int vtx) {
    for (int dof = 0; dof < 3; ++dof) m_scriptedGoal.erase(dof + vtx * 4);
  }

  void ungoalTheta(int vtx) { m_scriptedGoal.erase(3 + vtx * 4); }

  void unfreezeVertices(int vtx) {
    for (int dof = 0; dof < 3; ++dof)
      m_scriptedDegreesOfFreedom.erase(dof + vtx * 4);
  }

  void unfreezeTheta(int vtx) { m_scriptedDegreesOfFreedom.erase(3 + vtx * 4); }

  void freezeVertices(int vtx) {
    for (int dof = 0; dof < 3; ++dof)
      m_scriptedDegreesOfFreedom[dof + vtx * 4] = 0.0;
  }

  void freezeTheta(int vtx) { m_scriptedDegreesOfFreedom[vtx * 4 + 3] = 0.0; }

  void freezeVerticesWithTwist(int vtx) {
    for (int dof = 0; dof < 4; ++dof)
      m_scriptedDegreesOfFreedom[dof + vtx * 4] = 0.0;
  }

  void setVertexDisplacement(int vtx, const Vec3x& vel) {
    m_scriptedDegreesOfFreedom[4 * vtx + 0] = vel[0];
    m_scriptedDegreesOfFreedom[4 * vtx + 1] = vel[1];
    m_scriptedDegreesOfFreedom[4 * vtx + 2] = vel[2];
  }

  void setThetaDisplacement(int vtx, Scalar vel) {
    m_scriptedDegreesOfFreedom[4 * vtx + 3] = vel;
  }

  void setVertexGoal(int vtx, const Vec3x& pos) {
    m_scriptedGoal[4 * vtx + 0] = pos[0];
    m_scriptedGoal[4 * vtx + 1] = pos[1];
    m_scriptedGoal[4 * vtx + 2] = pos[2];
  }

  void setThetaGoal(int vtx, Scalar pos) { m_scriptedGoal[4 * vtx + 3] = pos; }

  Vec3x getVertexGoal(int vtx) const {
    Vec3x v = Vec3x::Zero();
    for (int i = 0; i < 3; ++i) {
      auto itr = m_scriptedGoal.find(4 * vtx + i);
      if (itr != m_scriptedGoal.end()) {
        v(i) = itr->second;
      }
    }

    return v;
  }

  Scalar getThetaGoal(int vtx) const {
    Scalar t = 0.0;
    auto itr = m_scriptedGoal.find(4 * vtx + 3);
    if (itr != m_scriptedGoal.end()) {
      t = itr->second;
    }
    return t;
  }

  void fixLHS(JacobianMatrixType& LHS) const;
  void fixRHS(VecXx& rhs) const;

  void fixRHS(JacobianMatrixType& LHS, VecXx& rhs, Scalar dt) const;

  void fixLHSAndRHS(JacobianMatrixType& LHS, VecXx& rhs, Scalar dt) const;

  void enforceDisplacements(VecXx& displacements, Scalar fraction = 1.0) const;
  void enforceVelocities(VecXx& velocities, Scalar dt) const;

  void makeBoundaryConditions(std::vector<VecXd>& bndry_normal, int ndof,
                              int offset = 0) const {
    for (auto dof = m_scriptedDegreesOfFreedom.begin();
         dof != m_scriptedDegreesOfFreedom.end(); ++dof) {
      bndry_normal.push_back(VecXd(ndof));
      bndry_normal.back()(dof->first + offset) = 1.0;
    }
  }

  size_t numScriptedDOF() const {
    return m_scriptedDegreesOfFreedom.size() + m_scriptedGoal.size();
  }

  void setRootFrame(const Vec3x& rootFrame) {
    m_rootFrame = rootFrame;
    m_enforceRootFrame = true;
  }

  void update(StrandState& geometry, const unsigned subStepId,
              const unsigned numberOfSubsteps);

  /**
   * @brief Computes rigid motion based on (scripted) root vertices
   * displacement.
   *
   * This could to be used either as initial guess in a Newton method, or as
   * safeguard motion if the regular solve fails. This DOFScriptingController
   * generates displacements: if the first edge is purely translated the
   * translation is propagated to all DOFs; otherwise a rigid motion based on
   * parallel transport of the first edge is applied.
   *
   * @param futureDOFs: position after rigid motion
   * @param currentDOFs: position before rigid motion
   */
  void computeRigidBodyMotion(VecXx& futureDOFs, const VecXx& currentDOFs);

  std::map<int, Scalar> m_scriptedDegreesOfFreedom;

  std::map<int, Scalar> m_scriptedGoal;

  bool m_enforceRootFrame;
  Vec3x m_rootFrame;
};

} /* namespace strandsim */
#endif /* VERTEXSCRIPTINGCONTROLLER_HH_ */
