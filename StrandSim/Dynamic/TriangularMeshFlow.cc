/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "TriangularMeshFlow.hh"

#include "../Collision/TriangularMesh.hh"
#include "../Dynamic/FluidScriptingController.hh"
#include "../Forces/GravitationForce.hh"
#include "../Utils/MathUtilities.hh"
#include "../Utils/ThreadUtils.hh"

namespace strandsim {
TriangularMeshFlow::TriangularMeshFlow(
    const std::shared_ptr<TriangularMesh>& mesh,
    const std::shared_ptr<FluidScriptingController>& controller,
    const Scalar& init_height)
    : m_mesh(mesh), m_controller(controller) {
  m_current_flow_height.resize(m_mesh->nv());
  m_future_flow_height.resize(m_mesh->nv());
  m_current_vertex_flow_velocity.resize(m_mesh->nv() * 3);
  m_touched.resize(m_mesh->nv());

  m_current_flow_velocity.resize(m_mesh->nf() * 3);
  m_future_flow_velocity.resize(m_mesh->nf() * 3);
  m_flow_mass_face.resize(m_mesh->nf());
  m_flow_rhs.resize(m_mesh->nf() * 3);
  m_flow_lhs.resize(m_mesh->nf() * 3);
  m_flow_cauchy_strain.resize(m_mesh->nf() * 3, 3);

  m_current_flow_velocity.setZero();
  m_future_flow_velocity.setZero();
  m_current_vertex_flow_velocity.setZero();

  const int num_components = m_controller->getNumComponents();

  m_current_flow_components.resize(m_mesh->nv() * num_components);
  m_current_flow_components.setZero();
  for (int i = 0; i < m_mesh->nv(); ++i) {
    m_current_flow_components(i * num_components) = 1.0;
  }
  m_future_flow_components = m_current_flow_components;

  //        for(int i = 0; i < m_mesh->nv(); ++i) {
  //            if(m_mesh->m_vertices[i](1) > 1.0) {
  //                m_current_flow_height(i) = init_height;
  //            } else {
  //                m_current_flow_height(i) = 0.0;
  //            }
  //        }
  //
  m_current_flow_height.setConstant(init_height);
  m_future_flow_height = m_current_flow_height;
  m_flow_mass_face.setZero();
  m_flow_rhs.setZero();
  m_flow_lhs.setOnes();

  for_each(0, (int)m_mesh->nf(), [this](int i) {
    m_flow_cauchy_strain.block<3, 3>(i * 3, 0).setIdentity();
  });

  m_shear_force.resize(m_mesh->nf() * 3);
  m_shear_force.setZero();

  const Scalar b = m_controller->getMeshFlowSlipLength();

  for (int i = 0; i < (int)m_mesh->nv(); ++i) {
    m_touched(i) = m_current_flow_height(i) >= b;
  }

  m_passing_vertices.resize(m_mesh->nv());
  m_passing_vertices.setZero();
}

const VecXx& TriangularMeshFlow::getCurrentFlowVelocity() const {
  return m_current_flow_velocity;
}

VecXx& TriangularMeshFlow::getCurrentFlowVelocity() {
  return m_current_flow_velocity;
}

const VecXx& TriangularMeshFlow::getCurrentFlowHeight() const {
  return m_current_flow_height;
}

VecXx& TriangularMeshFlow::getCurrentFlowHeight() {
  return m_current_flow_height;
}

void TriangularMeshFlow::check_passing_vertices() {
  if (!m_controller) return;

  const int num_verts = m_mesh->nv();
  const Scalar dx = m_controller->getDx();

  m_passing_vertices.setZero();

  VecXx distances(num_verts);

  for_each(0, num_verts, [&](int i) {
    const Vec3x pos = m_mesh->getVertex(i);
    distances[i] = m_controller->get_dense_liquid_signed_distance(pos);
  });

  for_each(0, (int)m_promising_vertices.size(), [&](int i) {
    int vert_idx = m_promising_vertices[i];
    m_passing_vertices[vert_idx] = distances[vert_idx] >= -1.5 * dx;
  });

  m_promising_vertices.reserve(num_verts);
  m_promising_vertices.resize(0);

  for (int i = 0; i < num_verts; ++i) {
    if (distances[i] < -0.5 * dx) {
      m_promising_vertices.push_back(i);
    }
  }
}

Vec3x TriangularMeshFlow::getFaceVelocity(int face_idx) const {
  return m_current_flow_velocity.segment<3>(face_idx * 3);
}

void TriangularMeshFlow::setFaceVelocity(int face_idx, const Vec3x& v) {
  m_current_flow_velocity.segment<3>(face_idx * 3) = v;
}

void TriangularMeshFlow::setVertexHeight(int local_idx, Scalar height) {
  m_current_flow_height(local_idx) = height;
}

Vec3x TriangularMeshFlow::getVertexVelocity(int local_idx) const {
  return m_current_vertex_flow_velocity.segment<3>(local_idx * 3);
}

Scalar TriangularMeshFlow::getVertexHeight(int local_idx) const {
  return m_current_flow_height(local_idx);
}

Scalar TriangularMeshFlow::getFaceHeight(int face_idx) const {
  const TriangularFace& f = m_mesh->getFace(face_idx);
  return (m_current_flow_height(f(0)) + m_current_flow_height(f(1)) +
          m_current_flow_height(f(2))) /
         3.0;
}

bool TriangularMeshFlow::isVertexPassing(int local_idx) const {
  return m_passing_vertices[local_idx];
}

void TriangularMeshFlow::step_dynamics(const Scalar dt) {
  Scalar total_vol = m_current_flow_height.sum();
  if (total_vol < 1e-8) return;

  advect_flow(dt);
  add_force(dt);
  update_flow_height(dt);

  {
    LockGuard lock(m_mesh->getGeometryMutex());
    m_current_flow_height = m_future_flow_height;
    m_current_flow_velocity = m_future_flow_velocity;
    m_current_flow_components = m_future_flow_components;
  }

  m_mesh->mapFaceVecToVertex(m_current_flow_velocity,
                             m_current_vertex_flow_velocity, 3);
}

void TriangularMeshFlow::setVertexComponents(int local_idx, const VecXx& c) {
  m_current_flow_components.segment(local_idx * c.size(), c.size()) = c;
}

const VecXx& TriangularMeshFlow::getCurrentComponents() const {
  return m_current_flow_components;
}

void TriangularMeshFlow::advect_flow(const Scalar dt) {
  // map face vel to vertices
  m_mesh->mapFaceVecToVertex(m_current_flow_velocity,
                             m_current_vertex_flow_velocity, 3);

  // advect velocity
  const int num_faces = m_mesh->nf();
  for_each(0, num_faces, [&](int i) {
    int target;
    Vec3x coord;
    m_mesh->getClosestBarycentric(
        i, -m_current_flow_velocity.segment<3>(i * 3) * dt, target, coord);
    const Vec3x& norm_target = m_mesh->m_face_frames.block<3, 1>(target * 3, 1);
    const Vec3x& norm_current = m_mesh->m_face_frames.block<3, 1>(i * 3, 1);

    const TriangularFace& f = m_mesh->m_faces[target];
    // sample from target
    Vec3x sampled_vel =
        m_current_vertex_flow_velocity.segment<3>(f(0) * 3) * coord(0) +
        m_current_vertex_flow_velocity.segment<3>(f(1) * 3) * coord(1) +
        m_current_vertex_flow_velocity.segment<3>(f(2) * 3) * coord(2);
    m_future_flow_velocity.segment<3>(i * 3) =
        parallelTransport(sampled_vel, norm_target, norm_current);
  });

  // advect height field
  const int num_verts = m_mesh->nv();
  const int num_components = m_controller->getNumComponents();
  for_each(0, num_verts, [&](int i) {
    const Vec3x& vel = m_current_vertex_flow_velocity.segment<3>(i * 3);
    int start_face_idx = m_mesh->getNeighborFaceAtDir(i, -vel * dt);

    const TriangularFace& f_start = m_mesh->m_faces[start_face_idx];
    const Vec3x offset_start =
        (m_mesh->m_vertices[f_start(0)] + m_mesh->m_vertices[f_start(1)] +
         m_mesh->m_vertices[f_start(2)]) /
            3.0 -
        m_mesh->m_vertices[i];

    int target;
    Vec3x coord;
    m_mesh->getClosestBarycentric(start_face_idx, -vel * dt - offset_start,
                                  target, coord);

    const TriangularFace& f_target = m_mesh->m_faces[target];

    Scalar h = m_current_flow_height(f_target(0)) * coord(0) +
               m_current_flow_height(f_target(1)) * coord(1) +
               m_current_flow_height(f_target(2)) * coord(2);
    m_future_flow_height(i) = h;

    VecXx c = m_current_flow_components.segment(f_target(0) * num_components,
                                                num_components) *
                  m_current_flow_height(f_target(0)) * coord(0) +
              m_current_flow_components.segment(f_target(1) * num_components,
                                                num_components) *
                  m_current_flow_height(f_target(1)) * coord(1) +
              m_current_flow_components.segment(f_target(2) * num_components,
                                                num_components) *
                  m_current_flow_height(f_target(2)) * coord(2);

    make_gibbs_simplex(c);

    m_future_flow_components.segment(i * num_components, num_components) = c;
  });
}

Scalar TriangularMeshFlow::getFaceVol(int face_idx) const {
  const TriangularFace& f = m_mesh->m_faces[face_idx];

  Scalar avg_fh = (m_current_flow_height(f(0)) + m_current_flow_height(f(1)) +
                   m_current_flow_height(f(2))) /
                  3.0;
  Scalar vol = m_mesh->m_face_areas(face_idx) * avg_fh;
  return vol;
}

Scalar TriangularMeshFlow::getVertexVol(int local_idx) const {
  return m_mesh->m_vertex_areas(local_idx) * m_current_flow_height(local_idx);
}

void TriangularMeshFlow::matrix_from_vec(const VecXx& vec,
                                         Eigen::SparseMatrix<Scalar>& mat) {
  mat.resize(vec.size(), vec.size());

  std::vector<Eigen::Triplet<Scalar> > tri_mat;
  tri_mat.reserve(vec.size());
  for (int i = 0; i < (int)vec.size(); i++) {
    if (vec(i) != 0.0) tri_mat.push_back(Eigen::Triplet<Scalar>(i, i, vec(i)));
  }
  mat.setFromTriplets(tri_mat.begin(), tri_mat.end());
}

void TriangularMeshFlow::matrix_from_vec3s(const VecXx& vec,
                                           Eigen::SparseMatrix<Scalar>& mat) {
  int num_terms = vec.size() / 3;
  mat.resize(num_terms, num_terms * 3);

  std::vector<Eigen::Triplet<Scalar> > tri_mat;
  tri_mat.reserve(vec.size());

  for (int i = 0; i < num_terms; ++i) {
    for (int r = 0; r < 3; ++r) {
      if (vec(i * 3 + r) != 0.0)
        tri_mat.push_back(Eigen::Triplet<Scalar>(i, i * 3 + r, vec(i * 3 + r)));
    }
  }

  mat.setFromTriplets(tri_mat.begin(), tri_mat.end());
}

void TriangularMeshFlow::update_flow_height(const Scalar dt) {
  const int num_faces = m_mesh->nf();
  const int num_verts = m_mesh->nv();

  VecXx div_u = m_mesh->m_div_operator * m_future_flow_velocity;

  for_each(0, num_verts,
           [&](int i) { m_future_flow_height(i) *= exp(-div_u(i) * dt); });

  if (m_controller->getMeshFlowEpsilon() > 0.0) {
    Eigen::SparseMatrix<Scalar> IL(num_verts, num_verts);
    IL.setIdentity();
    IL += m_mesh->m_G_v * m_mesh->m_laplacian_operator *
          m_controller->getMeshFlowEpsilon();

    Eigen::ConjugateGradient<Eigen::SparseMatrix<Scalar> > cg;
    cg.compute(IL);
    cg.setTolerance(1e-12);
    cg.setMaxIterations(1000);

    m_future_flow_height =
        cg.solveWithGuess(m_future_flow_height, m_future_flow_height);
  }

  MatXx vert_norms(num_verts, 3);
  for_each(0, num_verts, [&](int i) {
    vert_norms.block<1, 3>(i, 0) = m_mesh->m_vertex_normals[i].transpose();
  });

  MatXx grad_norm = m_mesh->m_grad_operator * vert_norms;

  const Scalar b = m_controller->getMeshFlowSlipLength();

  VecXx total_vol_vec(num_verts);
  VecXx new_total_vol_vec(num_verts);

  const Scalar max_flow_h = m_controller->getMeshFlowMaxHeight();

  for_each(0, num_verts, [&](int i) {
    Scalar bh = m_touched[i] ? b : 0.0;
    Scalar h = std::max(0.0, m_current_flow_height(i) - bh);
    Scalar vol = m_mesh->m_vertex_areas(i) * h;
    total_vol_vec(i) = vol;

    Scalar future_h = std::max(0.0, m_future_flow_height(i) - bh);
    Scalar future_vol = m_mesh->m_vertex_areas(i) * future_h;
    new_total_vol_vec(i) = future_vol;
  });

  Scalar new_total_vol = new_total_vol_vec.sum();
  Scalar total_vol = total_vol_vec.sum();

  Scalar prop = (new_total_vol == 0.0) ? 1.0 : (total_vol / new_total_vol);

  for_each(0, num_verts, [&](int i) {
    Scalar bh = m_touched[i] ? b : 0.0;
    m_future_flow_height(i) =
        std::min(max_flow_h,
                 std::max(0.0, (m_future_flow_height(i) - bh) * prop)) +
        bh;
    m_touched(i) = m_future_flow_height(i) >= b;
  });

  VecXx vertex_flow_vel(m_mesh->nv() * 3);

  m_mesh->mapFaceVecToVertex(m_future_flow_velocity, vertex_flow_vel, 3);

  MatXx mat_vert_flow_vel(m_mesh->nv(), 3);
  for_each(0, num_verts, [&](int i) {
    mat_vert_flow_vel.block<1, 3>(i, 0) =
        vertex_flow_vel.segment<3>(i * 3).transpose();
  });

  MatXx mat_grad_face_flow_vel = m_mesh->m_grad_operator * mat_vert_flow_vel;

  const int num_components = m_controller->getNumComponents();

  for_each(0, num_faces, [&](int i) {
    const TriangularFace& f = m_mesh->m_faces[i];

    if (m_future_flow_height(f(0)) < 1e-20 &&
        m_future_flow_height(f(1)) < 1e-20 &&
        m_future_flow_height(f(2)) < 1e-20) {
      m_flow_cauchy_strain.block<3, 3>(i * 3, 0).setIdentity();
      return;
    }

    const Mat3x gradu =
        mat_grad_face_flow_vel.block<3, 3>(i * 3, 0).transpose();
    const Mat3x& be_bar = m_flow_cauchy_strain.block<3, 3>(i * 3, 0);

    const Mat3x ff = Mat3x::Identity() + gradu * dt;
    const Scalar det_f = ff.determinant();

    const Scalar f_bar_factor = pow(det_f, -1.0 / 3.0);
    const Mat3x f_bar = ff * f_bar_factor;

    Mat3x be_bar_trial = f_bar * be_bar * f_bar.transpose();

    const Scalar factor = 1.0 / pow(be_bar_trial.determinant(), 1.0 / 3.0);
    be_bar_trial *= factor;

    VecXx color =
        m_future_flow_components.segment(f(0) * num_components,
                                         num_components) +
        m_future_flow_components.segment(f(1) * num_components,
                                         num_components) +
        m_future_flow_components.segment(f(2) * num_components, num_components);

    make_gibbs_simplex(color);

    Vec3x proj_func = m_controller->computeProjFunc(
        be_bar_trial, m_controller->getYieldStress(color),
        m_controller->getShearModulus(color), m_controller->getViscosity(color),
        m_controller->getFlowBehaviorIndex(color));

    const Scalar Ie_bar = be_bar_trial.trace() / 3.0;

    const Mat3x new_be_bar =
        proj_func(0) * (be_bar_trial - Ie_bar * Mat3x::Identity()) +
        Ie_bar * Mat3x::Identity();

    m_flow_cauchy_strain.block<3, 3>(i * 3, 0) = new_be_bar;
  });
}

void TriangularMeshFlow::transform_base() {
  const int num_faces = m_mesh->nf();
  for_each(0, num_faces, [this](int i) {
    Vec3x center;
    Vec3x translate;
    Eigen::Quaternion<Scalar> rot;
    m_mesh->getRigidTransform(i, center, translate, rot);

    m_current_flow_velocity.segment<3>(i * 3) =
        rot * Vec3x(m_current_flow_velocity.segment<3>(i * 3));
    m_future_flow_velocity.segment<3>(i * 3) =
        rot * Vec3x(m_future_flow_velocity.segment<3>(i * 3));

    Mat3x mrot = rot.toRotationMatrix();
    m_flow_cauchy_strain.block<3, 3>(i * 3, 0) =
        mrot * m_flow_cauchy_strain.block<3, 3>(i * 3, 0) * mrot.transpose();
  });
}

void TriangularMeshFlow::solve_flow_velocity() {
  const int num_faces = m_mesh->nf();

  for_each(0, num_faces, [this](int i) {
    if (m_flow_lhs(i * 3) > 1e-12) {
      m_future_flow_velocity.segment<3>(i * 3) =
          Vec3x(m_flow_rhs.segment<3>(i * 3).array() /
                m_flow_lhs.segment<3>(i * 3).array());
    } else {
      m_future_flow_velocity.segment<3>(i * 3).setZero();
    }

    Vec3x normal = m_mesh->m_face_frames.block<3, 1>(i * 3, 1);
    m_future_flow_velocity.segment<3>(i * 3) =
        (Mat3x::Identity() - normal * normal.transpose()) *
        m_future_flow_velocity.segment<3>(i * 3);
  });
}

void TriangularMeshFlow::add_force(const Scalar dt) {
  // compute mass
  const int num_faces = m_mesh->nf();
  const int num_verts = m_mesh->nv();

  //        VecXx grad_height = m_mesh->m_grad_operator * m_future_flow_height;

  MatXx mu_dev_b(num_faces * 3, 3);

  const int num_components = m_controller->getNumComponents();

  for_each(0, num_faces, [&](int i) {
    const TriangularFace& f = m_mesh->m_faces[i];

    Scalar h_avg = (m_future_flow_height(f(0)) + m_future_flow_height(f(1)) +
                    m_future_flow_height(f(2))) /
                   3.0;

    VecXx color = m_current_flow_components.segment(f(0) * num_components,
                                                    num_components) +
                  m_current_flow_components.segment(f(1) * num_components,
                                                    num_components) +
                  m_current_flow_components.segment(f(2) * num_components,
                                                    num_components);

    make_gibbs_simplex(color);

    m_flow_mass_face(i) =
        m_mesh->m_face_areas(i) * h_avg * m_controller->getDensity(color);

    mu_dev_b.block<3, 3>(i * 3, 0) =
        (m_flow_cauchy_strain.block<3, 3>(i * 3, 0) -
         Mat3x::Identity() *
             m_flow_cauchy_strain.block<3, 3>(i * 3, 0).trace() / 3.0) *
        m_controller->getShearModulus(color);
  });

  MatXx div_sigma =
      m_mesh->m_weight_vertex_face * m_mesh->m_div_operator * mu_dev_b;

  const Scalar dx = m_controller->getDx();

  // build rhs
  for_each(0, num_faces, [&](int i) {
    const TriangularFace& f = m_mesh->m_faces[i];
    Vec3x cur_pos = (m_mesh->m_vertices[f(0)] + m_mesh->m_vertices[f(1)] +
                     m_mesh->m_vertices[f(2)]) /
                    3.0;
    Scalar dist = m_controller->get_dense_liquid_signed_distance(cur_pos);
    if (dist < -1.0 * dx) {
      Vec3x vel = m_controller->get_velocity(cur_pos);
      Vec3x normal = m_mesh->m_face_frames.block<3, 1>(i * 3, 1);

      m_flow_rhs.segment<3>(i * 3) =
          (Mat3x::Identity() - normal * normal.transpose()) * vel *
          m_flow_mass_face(i);
    } else {
      const TriangularFace& f = m_mesh->m_faces[i];
      Vec3x cur_disp =
          (m_mesh->m_displacements[f(0)] + m_mesh->m_displacements[f(1)] +
           m_mesh->m_displacements[f(2)]) /
          3.0;
      Vec3x prev_disp = (m_mesh->m_previous_displacements[f(0)] +
                         m_mesh->m_previous_displacements[f(1)] +
                         m_mesh->m_previous_displacements[f(2)]) /
                        3.0;

      Scalar h_avg = (m_future_flow_height(f(0)) + m_future_flow_height(f(1)) +
                      m_future_flow_height(f(2))) /
                     3.0;
      Scalar vol = m_mesh->m_face_areas(i) * h_avg;

      Vec3x normal = m_mesh->m_face_frames.block<3, 1>(i * 3, 1);
      Vec3x shear_force = div_sigma.block<1, 3>(i, 0).transpose();
      m_shear_force.segment<3>(i * 3) = shear_force;

      Vec3x gt =
          (GravitationForce::getGravity() * dt - (cur_disp - prev_disp) / dt) *
              m_flow_mass_face(i) +
          shear_force * vol * dt;

      //            Vec3x gh = grad_height.segment<3>( i * 3 );

      m_flow_rhs.segment<3>(i * 3) =
          (m_future_flow_velocity.segment<3>(i * 3) * m_flow_mass_face(i) +
           (Mat3x::Identity() - normal * normal.transpose()) * gt);
    }
  });

  const Scalar b = m_controller->getMeshFlowSlipLength();

  // build lhs
  for_each(0, num_faces, [&](int i) {
    const TriangularFace& f = m_mesh->m_faces[i];
    m_flow_lhs.segment<3>(i * 3).setConstant(m_flow_mass_face(i));
  });

  // solve flow for precondition
  solve_flow_velocity();

  for_each(0, num_faces, [&](int i) {
    const TriangularFace& f = m_mesh->m_faces[i];
    Scalar h_avg = (m_future_flow_height(f(0)) + m_future_flow_height(f(1)) +
                    m_future_flow_height(f(2))) /
                   3.0;
    Scalar vol = m_mesh->m_face_areas(i) * h_avg;
    Scalar coeff = m_mesh->m_face_areas(i) / (b + h_avg / 3.0);

    VecXx color = m_current_flow_components.segment(f(0) * num_components,
                                                    num_components) +
                  m_current_flow_components.segment(f(1) * num_components,
                                                    num_components) +
                  m_current_flow_components.segment(f(2) * num_components,
                                                    num_components);

    make_gibbs_simplex(color);

    const Scalar eta = m_controller->getViscosity(color);
    const Scalar tau_Y = m_controller->getYieldStress(color) * 0.8164965809;
    const Scalar power = 1.0 - m_controller->getFlowBehaviorIndex(color);

    Scalar hu =
        h_avg / std::max(m_future_flow_velocity.segment<3>(i * 3).norm(), 1e-7);
    if (power < 0.0) hu = std::max(hu, 1e-10);

    const Scalar tilde_eta = tau_Y * hu + eta * pow(hu, power);

    m_flow_lhs.segment<3>(i * 3) += Vec3x::Constant(dt * coeff * tilde_eta);
  });

  // solve flow again
  solve_flow_velocity();
}
} /* namespace strandsim */
