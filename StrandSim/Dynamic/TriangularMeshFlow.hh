/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TRIANGULARMESHFLOW_HH_
#define TRIANGULARMESHFLOW_HH_

#include "../Collision/TriangularMesh.hh"
#include "../Core/Definitions.hh"
#include "../Utils/ThreadUtils.hh"

namespace strandsim {
class FluidScriptingController;

class TriangularMeshFlow {
  std::shared_ptr<TriangularMesh> m_mesh;

  // variables on vertex
  VecXx m_current_flow_height;
  VecXx m_future_flow_height;

  VecXx m_current_flow_components;
  VecXx m_future_flow_components;

  // variables on face
  VecXx m_current_flow_velocity;
  VecXx m_current_vertex_flow_velocity;
  VecXx m_future_flow_velocity;
  MatXx m_flow_cauchy_strain;  // 3n x 3

  VecXx m_flow_mass_face;
  VecXx m_flow_rhs;

  VecXx m_flow_lhs;

  VecXuc m_touched;

  std::vector<int> m_promising_vertices;
  VecXuc m_passing_vertices;

  std::shared_ptr<FluidScriptingController> m_controller;

  void matrix_from_vec3s(const VecXx& vec, Eigen::SparseMatrix<Scalar>& mat);
  void matrix_from_vec(const VecXx& vec, Eigen::SparseMatrix<Scalar>& mat);

 public:
  VecXx m_shear_force;

  TriangularMeshFlow(
      const std::shared_ptr<TriangularMesh>& mesh,
      const std::shared_ptr<FluidScriptingController>& controller,
      const Scalar& init_height);
  void advect_flow(const Scalar dt);
  void update_flow_height(const Scalar dt);
  void add_force(const Scalar dt);
  void solve_flow_velocity();
  void step_dynamics(const Scalar dt);
  void transform_base();

  void check_passing_vertices();

  const VecXx& getCurrentFlowHeight() const;
  VecXx& getCurrentFlowHeight();
  const VecXx& getCurrentFlowVelocity() const;
  VecXx& getCurrentFlowVelocity();
  const VecXx& getCurrentComponents() const;

  bool isVertexPassing(int local_idx) const;
  Scalar getVertexHeight(int local_idx) const;
  Vec3x getVertexVelocity(int local_idx) const;
  Scalar getFaceHeight(int face_idx) const;
  Scalar getFaceVol(int face_idx) const;
  Scalar getVertexVol(int local_idx) const;

  Vec3x getFaceVelocity(int face_idx) const;
  void setFaceVelocity(int face_idx, const Vec3x& v);

  void setVertexHeight(int local_idx, Scalar height);
  void setVertexComponents(int local_idx, const VecXx& c);
};
} /* namespace strandsim */
#endif /* TRIANGULARMESHFLOW_HH_ */
