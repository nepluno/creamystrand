/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "DOFScriptingController.hh"

#include "../Core/BandMatrix.hh"
#include "../Core/ElasticStrandUtils.hh"
#include "../Core/StrandState.hh"
#include "../Utils/MathUtilities.hh"
#include "../Utils/TextLog.hh"

namespace strandsim {

DOFScriptingController::DOFScriptingController() : m_enforceRootFrame(false) {
  m_scriptedDegreesOfFreedom.clear();
}

DOFScriptingController::DOFScriptingController(const Vec3xArray&)
    : m_enforceRootFrame(false) {
  m_scriptedDegreesOfFreedom.clear();
}

DOFScriptingController::~DOFScriptingController() {}

void DOFScriptingController::fixLHS(JacobianMatrixType& LHS) const {
  for (std::map<int, Scalar>::const_iterator dof =
           m_scriptedDegreesOfFreedom.begin();
       dof != m_scriptedDegreesOfFreedom.end(); ++dof)
    LHS.fixFirstDOFs<1>(dof->first);
}

void DOFScriptingController::fixRHS(VecXx& rhs) const {
  for (std::map<int, Scalar>::const_iterator dof =
           m_scriptedDegreesOfFreedom.begin();
       dof != m_scriptedDegreesOfFreedom.end(); ++dof)
    rhs[dof->first] = 0;
}

void DOFScriptingController::fixRHS(JacobianMatrixType& LHS, VecXx& rhs,
                                    Scalar dt) const {
  VecXx velocities(VecXx::Zero(rhs.rows()));

  for (std::map<int, Scalar>::const_iterator dof =
           m_scriptedDegreesOfFreedom.begin();
       dof != m_scriptedDegreesOfFreedom.end(); ++dof) {
    velocities[dof->first] = dof->second / dt;
  }

  LHS.multiply(rhs, -1, velocities);

  for (std::map<int, Scalar>::const_iterator dof =
           m_scriptedDegreesOfFreedom.begin();
       dof != m_scriptedDegreesOfFreedom.end(); ++dof) {
    rhs[dof->first] = velocities[dof->first];
  }
}

void DOFScriptingController::enforceDisplacements(VecXx& displacements,
                                                  Scalar fraction) const {
  for (std::map<int, Scalar>::const_iterator dof =
           m_scriptedDegreesOfFreedom.begin();
       dof != m_scriptedDegreesOfFreedom.end(); ++dof) {
    // std::cout << "Displacing " << dof->first << " by " << dof->second <<
    // '\n';
    displacements[dof->first] = dof->second * fraction;
  }
}

void DOFScriptingController::enforceVelocities(VecXx& velocities,
                                               Scalar dt) const {
  for (std::map<int, Scalar>::const_iterator dof =
           m_scriptedDegreesOfFreedom.begin();
       dof != m_scriptedDegreesOfFreedom.end(); ++dof) {
    //        std::cout << "Displacing " << dof->first << " by " << dof->second
    //        << '\n';
    velocities[dof->first] = dof->second / dt;
  }
}

void DOFScriptingController::fixLHSAndRHS(JacobianMatrixType& LHS, VecXx& rhs,
                                          Scalar dt) const {
  VecXx velocities(VecXx::Zero(rhs.rows()));

  for (std::map<int, Scalar>::const_iterator dof =
           m_scriptedDegreesOfFreedom.begin();
       dof != m_scriptedDegreesOfFreedom.end(); ++dof) {
    velocities[dof->first] = dof->second / dt;
  }

  LHS.multiply(rhs, -1, velocities);

  fixLHS(LHS);

  for (std::map<int, Scalar>::const_iterator dof =
           m_scriptedDegreesOfFreedom.begin();
       dof != m_scriptedDegreesOfFreedom.end(); ++dof) {
    rhs[dof->first] = velocities[dof->first];
  }
}

void DOFScriptingController::update(StrandState& geometry,
                                    const unsigned subStepId,
                                    const unsigned numberOfSubsteps) {}

void DOFScriptingController::computeRigidBodyMotion(VecXx& futureDOFs,
                                                    const VecXx& currentDOFs) {
  futureDOFs.resize(currentDOFs.size());
  VecXx displacements(currentDOFs.size());
  enforceDisplacements(displacements);

  const int m_ndof = currentDOFs.size();
  const Vec3x p0 = currentDOFs.segment<3>(0);
  const Vec3x p1 = currentDOFs.segment<3>(4);
  const Vec3x w0 = displacements.segment<3>(0);
  const Vec3x w1 = displacements.segment<3>(4);
  const Vec3x q0 = p0 + w0;
  const Vec3x q1 = p1 + w1;

  if (!isSmall(square((q1 - q0).squaredNorm() - (p1 - p0).squaredNorm())))
    ErrorStream(g_log, "") << "First edge length is not constant, diff "
                           << fabs((q1 - q0).squaredNorm() -
                                   (p1 - p0).squaredNorm());

  Vec3x u = (p1 - p0);
  const Scalar un = u.norm();
  if (!isSmall(un)) u /= un;

  Vec3x v = (q1 - q0);
  const Scalar vn = v.norm();
  if (!isSmall(vn)) v /= vn;

  if (isSmall(u.cross(v).squaredNorm()))  // pure translation
  {
    for (int i = 0; i < m_ndof; i += 4) {
      futureDOFs.segment<3>(i) = currentDOFs.segment<3>(i) + w0;
    }
  } else  // rigid motion
  {
    for (int i = 0; i < m_ndof; i += 4) {
      const Vec3x pi = currentDOFs.segment<3>(i);
      futureDOFs.segment<3>(i) = q0 + parallelTransport(pi - p0, u, v);
    }
  }
  // Finally just copy the thetas
  for (int i = 7; i < m_ndof; i += 4) {
    futureDOFs[i] = currentDOFs[i];
  }
}

} /* namespace strandsim */
