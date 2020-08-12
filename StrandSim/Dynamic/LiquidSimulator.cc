/**
 * \copyright 2019 Yun (Raymond) Fei, 2015 Yonghao Yue
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "LiquidSimulator.hh"

#include <deque>
#include <fstream>

#include "../Collision/TriangularMesh.hh"
#include "../Core/ElasticStrand.hh"
#include "../Forces/GravitationForce.hh"
#include "../Utils/SpherePattern.hh"
#include "AlgebraicMultigrid.hh"
#include "GeometricLevelGen.hh"
#include "ImplicitStepper.hh"
#include "StrandDynamicTraits.hh"
#include "TriangularMeshFlow.hh"

namespace strandsim {
LiquidSimulator::LiquidSimulator(std::vector<DistanceFieldObject>& fields_,
                                 std::vector<DistanceFieldObject>& sources_,
                                 std::vector<DistanceFieldObject>& terminators_,
                                 const std::vector<ElasticStrand*>& strands_,
                                 Scalar bucket_size_, int num_cells_,
                                 const Scalar dt_, const LiquidInfo li_)
    : FluidScriptingController(0, dt_),
      m_bucket_size(bucket_size_),
      m_num_nodes(num_cells_),
      m_liquid_info(li_),
      m_fields(fields_),
      m_sources(sources_),
      m_terminators(terminators_),
      m_strands(strands_) {
  sphere_pattern::generateSpherePattern(m_sphere_pattern);
  m_shooting_vol_accum.assign(sources_.size(), 0.0);

  const int num_strands = (int)m_strands.size();

  if (num_strands > 0) {
    m_strands_local_global_base.resize(num_strands);

    vector<int> tmp_buf(num_strands);

    for_each(0, num_strands,
             [&](int i) { tmp_buf[i] = m_strands[i]->getNumVertices(); });

    std::partial_sum(tmp_buf.begin(), tmp_buf.end(), tmp_buf.begin());

    for_each(0, num_strands, [&](int i) {
      if (i == 0)
        m_strands_local_global_base[i] = 0;
      else
        m_strands_local_global_base[i] = tmp_buf[i - 1];
    });

    const int num_vtx = tmp_buf[tmp_buf.size() - 1];

    m_strands_global_local.resize(num_vtx);

    for_each(0, num_strands, [&](int i) {
      const int num_vtx_local = m_strands[i]->getNumVertices();
      for (int j = 0; j < num_vtx_local; ++j) {
        const int idx_vtx = m_strands_local_global_base[i] + j;
        m_strands_global_local[idx_vtx] = pair<int, int>(i, j);
      }
    });

    m_vertex_weights.resize(num_vtx);
    m_strands_drag_coeffs.resize(num_vtx);
    m_strands_drag_coeffs.setZero();

    m_vertex_nodes_x.resize(num_vtx);
    m_vertex_nodes_y.resize(num_vtx);
    m_vertex_nodes_z.resize(num_vtx);
    m_vertex_nodes_p.resize(num_vtx);

    m_vertex_weights.resize(num_vtx);

    m_elastic_liquid_diffs.resize(num_vtx * 3);
    m_elastic_liquid_diffs.setZero();

    m_strands_submerged.resize(num_strands);
    m_strands_submerged.assign(num_strands, 0U);
  }

  const int num_fields = (int)m_fields.size();

  if (num_fields > 0) {
    m_mesh_local_global_base.resize(num_fields);

    vector<int> tmp_buf(num_fields);

    for_each(0, num_fields, [&](int i) {
      assert(m_fields[i].mesh_controller);

      std::shared_ptr<TriangularMesh> mesh =
          m_fields[i].mesh_controller->getCurrentMesh();
      assert(mesh);

      tmp_buf[i] = mesh->nv();
    });

    std::partial_sum(tmp_buf.begin(), tmp_buf.end(), tmp_buf.begin());

    for_each(0, num_fields, [&](int i) {
      if (i == 0)
        m_mesh_local_global_base[i] = 0;
      else
        m_mesh_local_global_base[i] = tmp_buf[i - 1];
    });

    const int num_vtx = tmp_buf[tmp_buf.size() - 1];

    m_mesh_global_local.resize(num_vtx);

    for_each(0, num_fields, [&](int i) {
      std::shared_ptr<TriangularMesh> mesh =
          m_fields[i].mesh_controller->getCurrentMesh();
      const int num_vtx_local = mesh->nv();
      for (int j = 0; j < num_vtx_local; ++j) {
        const int idx_vtx = m_mesh_local_global_base[i] + j;
        m_mesh_global_local[idx_vtx] = pair<int, int>(i, j);
      }
    });
  }
}

LiquidSimulator::~LiquidSimulator() {}

bool LiquidSimulator::isStrandSubmerged(int strand_idx) {
  return m_strands_submerged[strand_idx];
}

bool LiquidSimulator::processParticles() {
  std::cout << "[sample liquid distance fields]" << std::endl;
  sampleLiquidDistanceFields(m_time + m_dt);

  std::cout << "[drip from reservoir]" << std::endl;
  dripFromReservoir();

  std::cout << "[terminate particles]" << std::endl;
  terminateParticles();

  std::cout << "[capture bulk liquid]" << std::endl;
  captureBulkLiquid();

  std::cout << "[transfer bulk liquid mesh]" << std::endl;
  transferBulkLiquidMesh();

  std::cout << "[transfer liquid flow grid]" << std::endl;
  transferLiquidFlowGrid();

  std::cout << "[update optimal volume]" << std::endl;
  updateOptiVolume();

  std::cout << "[split liquid particles]" << std::endl;
  splitLiquidParticles();

  std::cout << "[merge liquid particles]" << std::endl;
  mergeLiquidParticles();

  std::cout << "[particle check weakened]" << std::endl;
  particleCheckWeakened();

  std::cout << "[solve color diffusion]" << std::endl;
  solveColorDiffusion();

  m_strands_drag_coeffs.setZero();

  return true;
}

bool LiquidSimulator::prepareDataStructure() {
  updateParticleBoundingBox();
  rebucketizeParticles();
  resampleNodes();
  return true;
}

void LiquidSimulator::initialize() {
  sampleLiquidDistanceFields(0.0);
  updateParticleBoundingBox();
  rebucketizeParticles();
  resampleNodes();
  computeWeight();
  updateSolidPhi();
  computeSignedDistance();
  correctLiquidParticles();
  updateSolidWeights();
  updateLiquidPhi();
  map_p2g();
  terminateParticles();
  updateOptiVolume();
  splitLiquidParticles();
  mergeLiquidParticles();
}
Scalar LiquidSimulator::getDx() const {
  return m_bucket_size / (Scalar)m_num_nodes;
}
Scalar LiquidSimulator::getInverseDCoeff() const {
  return inverse_D_coeff(getDx(), 2);
}
Scalar LiquidSimulator::getDensityAtPos(const Vec3x& pos) const {
  const Scalar vol = get_volume(pos);
  if (vol < 1e-20)
    return 0.0;
  else
    return get_mass(pos) / vol;
}

VecXx LiquidSimulator::getComponentsAtPos(const Vec3x& pos) const {
  const Scalar dx = getDx();
  VecXx color = interpolateValue(pos, m_node_components,
                                 m_grid_mincorner + Vec3x(0.5, 0.5, 0.5) * dx,
                                 VecXx::Zero(m_liquid_info.num_components), 1,
                                 m_liquid_info.num_components);

  make_gibbs_simplex(color);
  return color;
}

Scalar LiquidSimulator::getDensity(const VecXx& color) const {
  return m_liquid_info.liquid_density.dot(color);
}

Scalar LiquidSimulator::getBulkModulus(const VecXx& color) const {
  return m_liquid_info.liquid_bulk_modulus.dot(color);
}

Scalar LiquidSimulator::getShearModulus(const VecXx& color) const {
  return m_liquid_info.liquid_shear_modulus.dot(color);
}
int LiquidSimulator::getMaxIters() const { return m_liquid_info.pcg_max_iters; }
Scalar LiquidSimulator::getCriterion() const {
  return m_liquid_info.pressure_pcg_criterion;
}
Scalar LiquidSimulator::getContactAngle() const {
  return m_liquid_info.rest_contact_angle;
}
Scalar LiquidSimulator::getSurfTension(const VecXx& color) const {
  return m_liquid_info.surf_tension_coeff.dot(color);
}
Scalar LiquidSimulator::getSolidYoungModulus() const {
  return m_liquid_info.levelset_young_modulus;
}
Scalar LiquidSimulator::getViscosity(const VecXx& color) const {
  return m_liquid_info.flow_consistency_index.dot(color);
}
Scalar LiquidSimulator::getFlowBehaviorIndex(const VecXx& color) const {
  return m_liquid_info.flow_behavior_index.dot(color);
}
Scalar LiquidSimulator::getYieldStress(const VecXx& color) const {
  return m_liquid_info.plastic_yield_stress.dot(color);
}
Scalar LiquidSimulator::getEpsilon(const Vec3x& pos) const {
  const Scalar dx = getDx();
  return 1.0 - interpolateValue(pos, m_node_elastic_vf_p,
                                m_grid_mincorner + Vec3x(0.5, 0.5, 0.5) * dx,
                                0.0);
}
Scalar LiquidSimulator::cfl() {
  const Scalar dx = getDx();

  Scalar max_vel = 1e-20;
  const int num_parts = numParticles();
  for (int i = 0; i < num_parts; ++i) {
    max_vel = std::max(max_vel, m_v.segment<3>(i * 3).squaredNorm());
  }

  return dx / sqrt(max_vel);
}
Vec3i LiquidSimulator::getNodeHandle(int node_idx) const {
  int iz = node_idx / (m_num_nodes * m_num_nodes);
  node_idx -= iz * m_num_nodes * m_num_nodes;
  int iy = node_idx / m_num_nodes;
  int ix = node_idx - iy * m_num_nodes;

  return Vec3i(ix, iy, iz);
}

int LiquidSimulator::getNodeIndex(const Vec3i& handle) const {
  return handle(2) * m_num_nodes * m_num_nodes + handle(1) * m_num_nodes +
         m_num_nodes;
}

Vec3x LiquidSimulator::getNodePosSolidPhi(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx])
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3);
  else
    return nodePosFromBucket(bucket_idx, node_idx, Vec3x::Zero());
}

Vec3x LiquidSimulator::nodePosFromBucket(int bucket_idx, int raw_node_idx,
                                         const Vec3x& offset) const {
  Vec3i handle = m_particle_buckets.bucket_handle(bucket_idx);
  Vec3x bucket_left_corner =
      m_grid_mincorner + Vec3x(handle(0) * m_bucket_size,
                               handle(1) * m_bucket_size,
                               handle(2) * m_bucket_size);
  int iz = raw_node_idx / (m_num_nodes * m_num_nodes);
  int ixy = raw_node_idx - iz * m_num_nodes * m_num_nodes;
  int iy = ixy / m_num_nodes;
  int ix = ixy - iy * m_num_nodes;
  Vec3x node_pos = bucket_left_corner + (Vec3x(ix, iy, iz) + offset) * getDx();

  return node_pos;
}

Vec3x LiquidSimulator::getNodePos(int bucket_idx, int node_idx, int r) const {
  const static Vec3x offsets[] = {Vec3x(0.0, 0.5, 0.5), Vec3x(0.5, 0.0, 0.5),
                                  Vec3x(0.5, 0.5, 0.0), Vec3x(0.5, 0.5, 0.5),
                                  Vec3x(0.0, 0.0, 0.0)};

  if (m_bucket_activated[bucket_idx]) {
    const Scalar dx = getDx();
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3) + offsets[r] * dx;
  } else {
    return nodePosFromBucket(bucket_idx, node_idx, offsets[r]);
  }
}

Vec3x LiquidSimulator::getNodePosX(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx]) {
    const Scalar dx = getDx();
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3) +
           Vec3x(0.0, 0.5, 0.5) * dx;
  } else {
    return nodePosFromBucket(bucket_idx, node_idx, Vec3x(0.0, 0.5, 0.5));
  }
}

Vec3x LiquidSimulator::getNodePosY(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx]) {
    const Scalar dx = getDx();
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3) +
           Vec3x(0.5, 0.0, 0.5) * dx;
  } else {
    return nodePosFromBucket(bucket_idx, node_idx, Vec3x(0.5, 0.0, 0.5));
  }
}

Vec3x LiquidSimulator::getNodePosZ(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx]) {
    const Scalar dx = getDx();
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3) +
           Vec3x(0.5, 0.5, 0.0) * dx;
  } else {
    return nodePosFromBucket(bucket_idx, node_idx, Vec3x(0.5, 0.5, 0.0));
  }
}

Vec3x LiquidSimulator::getNodePosP(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx]) {
    const Scalar dx = getDx();
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3) +
           Vec3x(0.5, 0.5, 0.5) * dx;
  } else {
    return nodePosFromBucket(bucket_idx, node_idx, Vec3x(0.5, 0.5, 0.5));
  }
}

int LiquidSimulator::getNumNodes(int bucket_idx) const {
  if (m_bucket_activated[bucket_idx])
    return m_num_nodes * m_num_nodes * m_num_nodes;
  else
    return 0;
}

Scalar LiquidSimulator::getMeshFlowEpsilon() const {
  return m_liquid_info.mesh_flow_epsilon;
}

Scalar LiquidSimulator::solveHerschelBulkley(Scalar s_tr, Scalar dt,
                                             Scalar inv_eta, Scalar inv_m,
                                             Scalar yieldStress,
                                             Scalar __mu__) {
  const Scalar yc = sqrt(2.0 / 3.0) * yieldStress;
  return bisection_root_finding(
      s_tr, yc, m_liquid_info.pressure_pcg_criterion,
      m_liquid_info.pcg_max_iters, [&](const Scalar& s) {
        Scalar inside = std::max(0.0, (s - yc) * inv_eta);
        Scalar inside_pwr = (inside == 0.0) ? 0.0 : pow(inside, inv_m);

        return s - s_tr + 2.0 * __mu__ * dt * inside_pwr;
      });
}

Scalar LiquidSimulator::getMeshFlowMaxHeight() const {
  return m_liquid_info.mesh_flow_max_height;
}

Scalar LiquidSimulator::getMeshFlowSlipLength() const {
  return m_liquid_info.mesh_flow_slip_length;
}

Scalar LiquidSimulator::get_solid_phi_gradient(const Vec3x& position,
                                               Vec3x& n) const {
  const Scalar dx = getDx();
  return interpolateValueAndGradient(n, position, m_node_solid_phi,
                                     m_grid_mincorner, 3.0 * dx);
}

Vec3x LiquidSimulator::get_solid_gradient(const Vec3x& position) const {
  const Scalar dx = getDx();
  return interpolateGradient(position, m_node_solid_phi, m_grid_mincorner,
                             3.0 * dx);
}

Mat3x LiquidSimulator::get_solid_hessian(const Vec3x& position) const {
  const Scalar dx = getDx();
  Mat3x m;
  m.row(0) = (get_solid_gradient(position + Vec3x(0.5, 0.0, 0.0) * dx) -
              get_solid_gradient(position - Vec3x(0.5, 0.0, 0.0) * dx)) /
             dx;
  m.row(1) = (get_solid_gradient(position + Vec3x(0.0, 0.5, 0.0) * dx) -
              get_solid_gradient(position - Vec3x(0.0, 0.5, 0.0) * dx)) /
             dx;
  m.row(2) = (get_solid_gradient(position + Vec3x(0.0, 0.0, 0.5) * dx) -
              get_solid_gradient(position - Vec3x(0.0, 0.0, 0.5) * dx)) /
             dx;

  m = (m + m.transpose()) * 0.5;

  return m;
}

Vec3x LiquidSimulator::get_pressure_gradient(const Vec3x& position) const {
  const Scalar dx = getDx();
  Vec3x v;
  v(0) = interpolateValue(position, m_node_pressure_grad_x,
                          m_grid_mincorner + Vec3x(0.0, 0.5, 0.5) * dx, 0.0);
  v(1) = interpolateValue(position, m_node_pressure_grad_y,
                          m_grid_mincorner + Vec3x(0.5, 0.0, 0.5) * dx, 0.0);
  v(2) = interpolateValue(position, m_node_pressure_grad_z,
                          m_grid_mincorner + Vec3x(0.5, 0.5, 0.0) * dx, 0.0);

  return v;
}

Mat3x LiquidSimulator::get_pressure_hessian(const Vec3x& position) const {
  const Scalar dx = getDx();
  Mat3x m;
  m.row(0) = (get_pressure_gradient(position + Vec3x(0.5, 0.0, 0.0) * dx) -
              get_pressure_gradient(position - Vec3x(0.5, 0.0, 0.0) * dx)) /
             dx;
  m.row(1) = (get_pressure_gradient(position + Vec3x(0.0, 0.5, 0.0) * dx) -
              get_pressure_gradient(position - Vec3x(0.0, 0.5, 0.0) * dx)) /
             dx;
  m.row(2) = (get_pressure_gradient(position + Vec3x(0.0, 0.0, 0.5) * dx) -
              get_pressure_gradient(position - Vec3x(0.0, 0.0, 0.5) * dx)) /
             dx;

  m = (m + m.transpose()) * 0.5;

  return m;
}

Vec3x LiquidSimulator::computeProjFunc(
    const Mat3x& be_bar_trial, const Scalar& plastic_yield_stress,
    const Scalar& liquid_shear_modulus, const Scalar& flow_consistency_index,
    const Scalar& flow_behavior_index) const {
  auto dgdsmu = [](const Scalar& s, const Scalar& sig, const Scalar& eta,
                   const Scalar& mu, const Scalar& n,
                   const Scalar& h) -> Vec3x {
    if (eta == 0.0) {
      return Vec3x(sig / s, -sig / (s * s), 0.0);
    } else if (fabs(n - 1.0) < 1e-2) {
      const Scalar t2 = 1.0 / s;
      const Scalar t3 = sig * t2;
      const Scalar t4 = 1.0 / (s * s);
      const Scalar t5 = 1.0 / eta;
      const Scalar t7 = h * mu * t5 * 2.0;
      const Scalar t6 = exp(-t7);
      const Scalar t8 = t3 - 1.0;

      return Vec3x(t3 - t6 * t8, -sig * t4 + sig * t4 * t6,
                   h * t5 * t6 * t8 * 2.0);
    } else {
      const Scalar t2 = 1.0 / n;
      const Scalar t3 = n - 1.0;
      const Scalar t4 = s - sig;
      const Scalar t5 = t2 * t3;
      const Scalar t6 = pow(t4, t5);
      const Scalar t7 = pow(eta, -t2);
      const Scalar t14 = h * mu * t2 * t3 * t7 * 2.0;
      const Scalar t8 = t6 - t14;
      const Scalar t9 = 1.0 / t3;
      const Scalar t10 = n * t9;
      const Scalar t11 = pow(t8, t10);
      const Scalar t12 = sig + t11;
      const Scalar t13 = 1.0 / s;
      const Scalar t15 = t10 - 1.0;
      const Scalar t16 = pow(t8, t15);

      return Vec3x(t12 * t13,
                   -1.0 / (s * s) * t12 + pow(t4, t5 - 1.0) * t13 * t16,
                   h * t7 * t13 * t16 * -2.0);
    }
  };

  const Scalar Ie_bar = be_bar_trial.trace() / 3.0;
  const Mat3x s_trial =
      liquid_shear_modulus * (be_bar_trial - Mat3x::Identity() * Ie_bar);

  const Scalar s_trial_len = s_trial.norm();
  const Scalar tilde_sigma_Y = sqrt(2.0 / 3.0) * plastic_yield_stress;
  const Scalar phi_trial = s_trial_len - tilde_sigma_Y;

  Mat3x new_be_bar;

  if (phi_trial <= 0.0) {
    return Vec3x(1.0, 0.0, 0.0);
  } else {
    const Scalar __mu__ = Ie_bar * liquid_shear_modulus;
    return dgdsmu(s_trial_len, tilde_sigma_Y, flow_consistency_index, __mu__,
                  flow_behavior_index, m_dt);
  }
}

void LiquidSimulator::acceptFutureStates() {
  m_b = m_b_plus;
  m_Fe = m_Fe_plus;
}

bool LiquidSimulator::map_g2p() {
  const int num_part = numParticles();

  const Scalar invD = getInverseDCoeff();

  for_each(0, num_part, [&](int pidx) {
    auto& indices_x = m_particle_nodes_x[pidx];
    auto& indices_y = m_particle_nodes_y[pidx];
    auto& indices_z = m_particle_nodes_z[pidx];

    auto& weights = m_particle_weights[pidx];
    const Vec3x& pos = m_x.segment<3>(pidx * 3);
    //
    //            const Mat27x3f& grad_x = m_particle_grad_x[pidx];
    //            const Mat27x3f& grad_y = m_particle_grad_y[pidx];
    //            const Mat27x3f& grad_z = m_particle_grad_z[pidx];

    m_v.segment<3>(pidx * 3).setZero();

    Vec3x fv = Vec3x::Zero();
    Vec3x gp = Vec3x::Zero();
    //            Mat3x B = Mat3x::Zero();
    Mat3x C = Mat3x::Zero();

    for (int i = 0; i < indices_x.rows(); ++i) {
      const int node_bucket_idx = indices_x(i, 0);
      if (!m_bucket_activated[node_bucket_idx]) continue;

      const int node_idx = indices_x(i, 1);

      Scalar fnv = m_node_vel_fluid_plus_x[node_bucket_idx](node_idx);

      Vec3x np = getNodePosX(node_bucket_idx, node_idx);

      fv(0) += fnv * weights(i, 0);

      gp(0) +=
          m_node_pressure_grad_x[node_bucket_idx](node_idx) * weights(i, 0);

      //                B.row(0) += fnv * grad_x.row(i);
      C.row(0) += fnv * weights(i, 0) * (np - pos).transpose() * invD;
    }

    for (int i = 0; i < indices_y.rows(); ++i) {
      const int node_bucket_idx = indices_y(i, 0);
      if (!m_bucket_activated[node_bucket_idx]) continue;

      const int node_idx = indices_y(i, 1);

      Scalar fnv = m_node_vel_fluid_plus_y[node_bucket_idx](node_idx);

      Vec3x np = getNodePosY(node_bucket_idx, node_idx);

      fv(1) += fnv * weights(i, 1);

      gp(1) +=
          m_node_pressure_grad_y[node_bucket_idx](node_idx) * weights(i, 1);
      //                B.row(1) += fnv * grad_y.row(i);
      C.row(1) += fnv * weights(i, 1) * (np - pos).transpose() * invD;
    }

    for (int i = 0; i < indices_z.rows(); ++i) {
      const int node_bucket_idx = indices_z(i, 0);
      if (!m_bucket_activated[node_bucket_idx]) continue;

      const int node_idx = indices_z(i, 1);

      Scalar fnv = m_node_vel_fluid_plus_z[node_bucket_idx](node_idx);

      Vec3x np = getNodePosZ(node_bucket_idx, node_idx);

      fv(2) += fnv * weights(i, 2);

      gp(2) +=
          m_node_pressure_grad_z[node_bucket_idx](node_idx) * weights(i, 2);

      //                B.row(2) += fnv * grad_z.row(i);
      C.row(2) += fnv * weights(i, 2) * (np - pos).transpose() * invD;
    }

    m_v.segment<3>(pidx * 3) = fv;

    const VecXx& color = m_components.segment(
        pidx * m_liquid_info.num_components, m_liquid_info.num_components);
    const Scalar liquid_shear_modulus = getShearModulus(color);

    if (liquid_shear_modulus > 1e-20) {
      const Mat3x f = Mat3x::Identity() + C * m_dt;
      // check_isnan("f", f.sum());

      const Scalar det_f = f.determinant();
      if (det_f < 1e-20 || std::isnan(det_f) || det_f > 1e+20) {
        // particle is invalid, kill it
        m_vol(pidx) = 0.0;
        return;
      }

      const Mat3x& be_bar = m_b.block<3, 3>(pidx * 3, 0);

      // check_isnan("be_bar", be_bar.sum());

      Mat3x F;
      Mat3x be_bar_trial;
      // Handle Tearing
      if (!m_weakened[pidx]) {
        F = f * m_Fe.block<3, 3>(pidx * 3, 0);
        const Scalar f_bar_factor = pow(det_f, -1.0 / 3.0);
        const Mat3x f_bar = f * f_bar_factor;

        be_bar_trial = f_bar * be_bar * f_bar.transpose();

        // check_isnan("F_0", F.sum());
        // check_isnan("be_bar_trial_0", be_bar_trial.sum());
      } else {
        // check the largest eigenvalue.

        const Scalar f_bar_factor = pow(det_f, -1.0 / 3.0);
        const Mat3x f_bar = f * f_bar_factor;

        // check_isnan("f_bar", f_bar.sum());

        be_bar_trial = f_bar * be_bar * f_bar.transpose();

        Eigen::SelfAdjointEigenSolver<Mat3x> es;
        es.compute(be_bar);
        const Vec3x Sp = es.eigenvalues();

        const Scalar largestPrevEigenValue =
            std::max(std::max(Sp(0), Sp(1)), Sp(2));

        es.compute(be_bar_trial);
        const Vec3x Sn = es.eigenvalues();
        const Scalar largestTrialEigenValue =
            std::max(std::max(Sn(0), Sn(1)), Sn(2));

        if (largestTrialEigenValue > largestPrevEigenValue) {
          // cancel the shear
          // compute polar decomposition of f_bar, and only use its rotation
          // part; i.e., compute svd of f_bar, and set the magnitudes of the
          // singular values to one (retaining its sign).

          Eigen::JacobiSVD<Mat3x> svd(
              f_bar, Eigen::ComputeFullU | Eigen::ComputeFullV);
          const Mat3x U = svd.matrixU();
          const Mat3x V = svd.matrixV();
          Mat3x S = svd.singularValues().asDiagonal();

          S(0, 0) = (S(0, 0) >= 0.0) ? 1.0 : -1.0;
          S(1, 1) = (S(1, 1) >= 0.0) ? 1.0 : -1.0;
          S(2, 2) = (S(2, 2) >= 0.0) ? 1.0 : -1.0;

          const Mat3x corrected_f_bar = U * S * V.transpose();
          const Mat3x corrected_f = corrected_f_bar / f_bar_factor;

          // check_isnan("corrected_f", corrected_f.sum());

          F = corrected_f * m_Fe.block<3, 3>(pidx * 3, 0);
          // check_isnan("F", F.sum());

          be_bar_trial = corrected_f_bar * be_bar * corrected_f_bar.transpose();
        } else {
          F = f * m_Fe.block<3, 3>(pidx * 3, 0);
        }

        const Scalar Jprev = m_J(pidx);
        const Scalar J = F.determinant();

        if (J < 1e-20 || std::isnan(J)) {
          // particle is invalid, kill it
          m_vol(pidx) = 0.0;
          return;
        }

        if ((J > Jprev) && (J > 1.0)) F *= pow(Jprev / J, 1.0 / 3.0);

        // check_isnan("F_1", F.sum());

        // check_isnan("be_bar_trial_1", be_bar_trial.sum());
      }

      const Scalar det_bbt = be_bar_trial.determinant();

      if (det_bbt < 1e-20 || std::isnan(det_bbt)) {
        // particle is invalid, kill it
        m_vol(pidx) = 0.0;
        return;
      }

      const Scalar factor = 1.0 / std::max(1e-4, pow(det_bbt, 1.0 / 3.0));
      be_bar_trial *= factor;

      // check_isnan("be_bar_trial_2", be_bar_trial.sum());

      // check_isnan("F_2", F.sum());

      // Apply Plasticity
      Vec3x proj_func = computeProjFunc(
          be_bar_trial, getYieldStress(color), liquid_shear_modulus,
          getViscosity(color), getFlowBehaviorIndex(color));

      // check_isnan("proj_func", proj_func.sum());

      const Scalar Ie_bar = be_bar_trial.trace() / 3.0;

      const Mat3x new_be_bar =
          proj_func(0) * (be_bar_trial - Ie_bar * Mat3x::Identity()) +
          Ie_bar * Mat3x::Identity();

      // check_isnan("new_be_bar", new_be_bar.sum());

      m_Fe_plus.block<3, 3>(pidx * 3, 0) = F;

      m_b_plus.block<3, 3>(pidx * 3, 0) = new_be_bar;
      m_b_trial.block<3, 3>(pidx * 3, 0) = be_bar_trial;
      m_proj_func.segment<3>(pidx * 3) = proj_func;

      // Update J
      m_J(pidx) = F.determinant();

      // kill over-stretched vertices
      if (m_J(pidx) < 1e-20) m_vol(pidx) = 0.0;
    } else {
      m_Fe_plus.block<3, 3>(pidx * 3, 0).setIdentity();
      m_b_plus.block<3, 3>(pidx * 3, 0).setIdentity();
      m_b_trial.block<3, 3>(pidx * 3, 0).setIdentity();
      m_proj_func.segment<3>(pidx * 3) = Vec3x(1, 0, 0);
      m_J(pidx) = 1.0;
    }

    // Damping
    const Scalar cstretch = 1.0 - m_liquid_info.affine_stretch_damping;
    const Scalar crotate = 1.0 - m_liquid_info.affine_rotate_damping;

    m_B.block<3, 3>(pidx * 3, 0) =
        ((cstretch + crotate) * C + (cstretch - crotate) * C.transpose()) * 0.5;

    m_v.segment<3>(pidx * 3) *= 1.0 - m_liquid_info.velocity_damping;
  });

  removeEmptyParticles();
  return true;
}

bool LiquidSimulator::updateDragCoeffNodes() {
  const Scalar dx = getDx();
  const Scalar dV = dx * dx * dx;

  m_vertex_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    m_node_drag_coeff_x[bucket_idx].setZero();
    m_node_drag_coeff_y[bucket_idx].setZero();
    m_node_drag_coeff_z[bucket_idx].setZero();

    const auto& bucket_node_vertex_x = m_node_vertex_x[bucket_idx];
    const auto& bucket_node_vertex_y = m_node_vertex_y[bucket_idx];
    const auto& bucket_node_vertex_z = m_node_vertex_z[bucket_idx];

    const int num_nodes = getNumNodes(bucket_idx);

    for (int i = 0; i < num_nodes; ++i) {
      // ignore air nodes
      if (m_node_mass_fluid_x[bucket_idx][i] < 1e-20) continue;

      const auto& node_vertex_x = bucket_node_vertex_x[i];

      Scalar D = 0.0;

      for (auto& pair : node_vertex_x) {
        const int pidx = pair.first;
        auto& weights = m_vertex_weights[pidx];

        const Scalar w = weights(pair.second, 0);

        D += w * m_strands_drag_coeffs(pidx);
      }

      m_node_drag_coeff_x[bucket_idx][i] = D;
    }

    for (int i = 0; i < num_nodes; ++i) {
      // ignore air nodes
      if (m_node_mass_fluid_y[bucket_idx][i] < 1e-20) continue;

      const auto& node_vertex_y = bucket_node_vertex_y[i];

      Scalar D = 0.0;

      for (auto& pair : node_vertex_y) {
        const int pidx = pair.first;
        auto& weights = m_vertex_weights[pidx];

        const Scalar w = weights(pair.second, 1);

        D += w * m_strands_drag_coeffs(pidx);
      }

      m_node_drag_coeff_y[bucket_idx][i] = D;
    }

    for (int i = 0; i < num_nodes; ++i) {
      // ignore air nodes
      if (m_node_mass_fluid_z[bucket_idx][i] < 1e-20) continue;

      const auto& node_vertex_z = bucket_node_vertex_z[i];

      Scalar D = 0.0;

      for (auto& pair : node_vertex_z) {
        const int pidx = pair.first;
        auto& weights = m_vertex_weights[pidx];

        const Scalar w = weights(pair.second, 2);

        D += w * m_strands_drag_coeffs(pidx);
      }

      m_node_drag_coeff_z[bucket_idx][i] = D;
    }
  });

  // check_isnan("drag_coeff", m_node_drag_coeff_x, m_node_drag_coeff_y,
  // m_node_drag_coeff_z);

  return true;
}

bool LiquidSimulator::map_p2g() {
  const Scalar dx = getDx();
  const Scalar dV = dx * dx * dx;
  const Scalar max_packing = 0.9068996821;  // Thue's Theorem on Circle Packing

  // check_isnan("x_before_p2g", m_x);
  // check_isnan("v_before_p2g", m_v);
  // check_isnan("m_before_p2g", m_m);

  m_b_plus = m_b;
  m_Fe_plus = m_Fe;

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    m_node_mass_fluid_x[bucket_idx].setZero();
    m_node_vel_fluid_x[bucket_idx].setZero();
    m_node_vol_fluid_x[bucket_idx].setZero();
    m_node_vel_elastic_x[bucket_idx].setZero();
    m_node_pressure_grad_x[bucket_idx].setZero();

    m_node_mass_fluid_y[bucket_idx].setZero();
    m_node_vel_fluid_y[bucket_idx].setZero();
    m_node_vol_fluid_y[bucket_idx].setZero();
    m_node_vel_elastic_y[bucket_idx].setZero();
    m_node_pressure_grad_y[bucket_idx].setZero();

    m_node_mass_fluid_z[bucket_idx].setZero();
    m_node_vel_fluid_z[bucket_idx].setZero();
    m_node_vol_fluid_z[bucket_idx].setZero();
    m_node_vel_elastic_z[bucket_idx].setZero();
    m_node_pressure_grad_z[bucket_idx].setZero();

    m_node_vol_change_p[bucket_idx].setOnes();
    m_node_components[bucket_idx].setZero();

    const auto& bucket_node_particles_x = m_node_particles_x[bucket_idx];
    const auto& bucket_node_particles_y = m_node_particles_y[bucket_idx];
    const auto& bucket_node_particles_z = m_node_particles_z[bucket_idx];
    const auto& bucket_node_particles_p = m_node_particles_p[bucket_idx];

    const auto& bucket_node_vertex_x = m_node_vertex_x[bucket_idx];
    const auto& bucket_node_vertex_y = m_node_vertex_y[bucket_idx];
    const auto& bucket_node_vertex_z = m_node_vertex_z[bucket_idx];
    const auto& bucket_node_vertex_p = m_node_vertex_p[bucket_idx];

    const int num_nodes = getNumNodes(bucket_idx);

    for (int i = 0; i < num_nodes; ++i) {
      const auto& node_particles_x = bucket_node_particles_x[i];
      const Vec3x& np = getNodePosX(bucket_idx, i);

      Scalar p_fluid = 0.0;
      Scalar mass_fluid = 0.0;
      Scalar vol_fluid = 0.0;
      Scalar gp_fluid = 0.0;
      Scalar bulk_fluid = 0.0;

      for (auto& pair : node_particles_x) {
        const int pidx = pair.first;

        auto& weights = m_particle_weights[pidx];
        const Vec3x& m = m_m.segment<3>(pidx * 3);
        const Vec3x& v = m_v.segment<3>(pidx * 3);
        const Vec3x& pos = m_x.segment<3>(pidx * 3);
        const Mat3x& B = m_B.block<3, 3>(pidx * 3, 0);
        const Scalar fvol = m_vol(pidx);

        const Scalar vel = v(0) + B.row(0).dot(np - pos);
        p_fluid += vel * m(0) * weights(pair.second, 0);
        mass_fluid += m(0) * weights(pair.second, 0);
        vol_fluid += fvol * weights(pair.second, 0);
      }

      if (mass_fluid > 1e-20) {
        m_node_vel_fluid_x[bucket_idx](i) = p_fluid / mass_fluid;
      }

      m_node_mass_fluid_x[bucket_idx](i) = mass_fluid;
      m_node_vol_fluid_x[bucket_idx](i) = vol_fluid;

      const auto& node_vertex_x = bucket_node_vertex_x[i];

      Scalar p_elastic = 0.0;
      Scalar m_elastic = 0.0;
      for (auto& pair : node_vertex_x) {
        const int pidx = pair.first;
        const int strand_idx = m_strands_global_local[pidx].first;
        const int local_idx = m_strands_global_local[pidx].second;
        auto& weights = m_vertex_weights[pidx];

        const Scalar w = weights(pair.second, 0);

        const Scalar mass = m_strands[strand_idx]->getVertexMass(local_idx);
        p_elastic += m_strands[strand_idx]->getStepper()->velocities()(
                         local_idx * 4 + 0) *
                     mass * w;
        m_elastic += mass * w;
      }

      if (m_elastic > 1e-20) {
        m_node_vel_elastic_x[bucket_idx](i) = p_elastic / m_elastic;
      }
    }

    for (int i = 0; i < num_nodes; ++i) {
      const auto& node_particles_y = bucket_node_particles_y[i];
      const Vec3x& np = getNodePosY(bucket_idx, i);

      Scalar p_fluid = 0.0;
      Scalar mass_fluid = 0.0;
      Scalar vol_fluid = 0.0;
      Scalar gp_fluid = 0.0;

      for (auto& pair : node_particles_y) {
        const int pidx = pair.first;

        auto& weights = m_particle_weights[pidx];
        const Vec3x& m = m_m.segment<3>(pidx * 3);
        const Vec3x& v = m_v.segment<3>(pidx * 3);
        const Vec3x& pos = m_x.segment<3>(pidx * 3);
        const Mat3x& B = m_B.block<3, 3>(pidx * 3, 0);
        const Scalar fvol = m_vol(pidx);

        const Scalar vel = v(1) + B.row(1).dot(np - pos);
        p_fluid += vel * m(1) * weights(pair.second, 1);
        mass_fluid += m(1) * weights(pair.second, 1);
        vol_fluid += fvol * weights(pair.second, 1);
      }

      if (mass_fluid > 1e-20) {
        m_node_vel_fluid_y[bucket_idx](i) = p_fluid / mass_fluid;
      }

      m_node_mass_fluid_y[bucket_idx](i) = mass_fluid;
      m_node_vol_fluid_y[bucket_idx](i) = vol_fluid;

      const auto& node_vertex_y = bucket_node_vertex_y[i];

      Scalar p_elastic = 0.0;
      Scalar m_elastic = 0.0;
      for (auto& pair : node_vertex_y) {
        const int pidx = pair.first;
        const int strand_idx = m_strands_global_local[pidx].first;
        const int local_idx = m_strands_global_local[pidx].second;
        auto& weights = m_vertex_weights[pidx];

        const Scalar w = weights(pair.second, 1);

        const Scalar mass = m_strands[strand_idx]->getVertexMass(local_idx);
        p_elastic += m_strands[strand_idx]->getStepper()->velocities()(
                         local_idx * 4 + 1) *
                     mass * w;
        m_elastic += mass * w;
      }

      if (m_elastic > 1e-20) {
        m_node_vel_elastic_y[bucket_idx](i) = p_elastic / m_elastic;
      }
    }

    for (int i = 0; i < num_nodes; ++i) {
      const auto& node_particles_z = bucket_node_particles_z[i];
      const Vec3x& np = getNodePosZ(bucket_idx, i);

      Scalar p_fluid = 0.0;
      Scalar mass_fluid = 0.0;
      Scalar vol_fluid = 0.0;
      Scalar gp_fluid = 0.0;

      for (auto& pair : node_particles_z) {
        const int pidx = pair.first;

        auto& weights = m_particle_weights[pidx];
        const Vec3x& m = m_m.segment<3>(pidx * 3);
        const Vec3x& v = m_v.segment<3>(pidx * 3);
        const Vec3x& pos = m_x.segment<3>(pidx * 3);
        const Mat3x& B = m_B.block<3, 3>(pidx * 3, 0);
        const Scalar fvol = m_vol(pidx);

        const Scalar vel = v(2) + B.row(2).dot(np - pos);
        p_fluid += vel * m(2) * weights(pair.second, 2);
        mass_fluid += m(2) * weights(pair.second, 2);
        vol_fluid += fvol * weights(pair.second, 2);
      }

      if (mass_fluid > 1e-20) {
        m_node_vel_fluid_z[bucket_idx](i) = p_fluid / mass_fluid;
      }

      m_node_mass_fluid_z[bucket_idx](i) = mass_fluid;
      m_node_vol_fluid_z[bucket_idx](i) = vol_fluid;

      const auto& node_vertex_z = bucket_node_vertex_z[i];

      Scalar p_elastic = 0.0;
      Scalar m_elastic = 0.0;
      for (auto& pair : node_vertex_z) {
        const int pidx = pair.first;
        const int strand_idx = m_strands_global_local[pidx].first;
        const int local_idx = m_strands_global_local[pidx].second;
        auto& weights = m_vertex_weights[pidx];

        const Scalar w = weights(pair.second, 2);

        const Scalar mass = m_strands[strand_idx]->getVertexMass(local_idx);
        p_elastic += m_strands[strand_idx]->getStepper()->velocities()(
                         local_idx * 4 + 2) *
                     mass * w;
        m_elastic += mass * w;
      }

      if (m_elastic > 1e-20) {
        m_node_vel_elastic_z[bucket_idx](i) = p_elastic / m_elastic;
      }
    }

    for (int i = 0; i < num_nodes; ++i) {
      const auto& node_particles_p = bucket_node_particles_p[i];
      const Vec3x& np = getNodePosP(bucket_idx, i);

      Scalar vol_change_fluid = 0.0;
      Scalar weight_fluid = 0.0;
      VecXx color_fluid = VecXx::Zero(m_liquid_info.num_components);

      for (auto& pair : node_particles_p) {
        const int pidx = pair.first;

        auto& weights = m_particle_weights[pidx];
        const Scalar J = clamp(m_J(pidx), 0.5, 2.0);
        const VecXx& color = m_components.segment(
            pidx * m_liquid_info.num_components, m_liquid_info.num_components);

        vol_change_fluid += J * weights(pair.second, 3);
        color_fluid += color * weights(pair.second, 3);
        weight_fluid += weights(pair.second, 3);
      }

      Scalar center_J;
      if (weight_fluid > 1e-20) {
        center_J = vol_change_fluid / weight_fluid;
      } else {
        center_J = 1.0;
      }
      m_node_vol_change_p[bucket_idx](i) = center_J;

      make_gibbs_simplex(color_fluid);

      m_node_components[bucket_idx].segment(i * m_liquid_info.num_components,
                                            m_liquid_info.num_components) =
          color_fluid;

      const auto& node_vertex_p = bucket_node_vertex_p[i];

      Scalar vol_elastic = 0.0;

      for (auto& pair : node_vertex_p) {
        const int pidx = pair.first;
        const int strand_idx = m_strands_global_local[pidx].first;
        const int local_idx = m_strands_global_local[pidx].second;

        auto& weights = m_vertex_weights[pidx];
        const Scalar len = m_strands[strand_idx]->getVoronoiLength(local_idx);
        const Scalar ra = m_strands[strand_idx]->getRadiusA(local_idx);
        const Scalar rb = m_strands[strand_idx]->getRadiusB(local_idx);
        const Scalar vol = M_PI * ra * rb * len;

        vol_elastic += vol * weights(pair.second, 3);
      }

      m_node_elastic_vf_p[bucket_idx][i] =
          clamp(vol_elastic / dV, 0.0, max_packing);
    }

    m_node_vel_fluid_plus_x[bucket_idx] = m_node_vel_fluid_x[bucket_idx];
    m_node_vel_fluid_plus_y[bucket_idx] = m_node_vel_fluid_y[bucket_idx];
    m_node_vel_fluid_plus_z[bucket_idx] = m_node_vel_fluid_z[bucket_idx];
  });

  // check_isnan("fluidvel_after_p2g", m_node_vel_fluid_plus_x,
  // m_node_vel_fluid_plus_y, m_node_vel_fluid_plus_z);
  // check_isnan("elastic_vf_p_after_p2g", m_node_elastic_vf_p);
  // check_isnan("vol_change_after_p2g", m_node_vol_change_p);

  return true;
}

void LiquidSimulator::constrainSurfFlowVelocity() {
  // gather back to interfacing and inside edges
  const int num_edges = m_strands_inside_segs.size();
  if (!num_edges) return;

  // check_isnan("flowvel_before_constraint", (int) m_strands.size(), [this]
  // (int i) {return m_strands[i]->getStepper()->getFutureFlowVelocity();});

  for_each(0, num_edges, [&](int eidx) {
    const int strand_idx = m_strands_inside_segs[eidx].first;
    const int local_idx = m_strands_inside_segs[eidx].second;

    auto& indices_x = m_inside_edge_nodes_x[eidx];
    auto& indices_y = m_inside_edge_nodes_y[eidx];
    auto& indices_z = m_inside_edge_nodes_z[eidx];

    auto& weights = m_inside_edge_weights[eidx];

    Vec3x fv = Vec3x::Zero();

    for (int i = 0; i < indices_x.rows(); ++i) {
      const int node_bucket_idx = indices_x(i, 0);
      if (node_bucket_idx < 0 || !m_bucket_activated[node_bucket_idx]) continue;

      const int node_idx = indices_x(i, 1);
      if (node_idx < 0) continue;

      Scalar fnv = m_node_vel_fluid_plus_x[node_bucket_idx](node_idx);

      fv(0) += fnv * weights(i, 0);
    }

    for (int i = 0; i < indices_y.rows(); ++i) {
      const int node_bucket_idx = indices_y(i, 0);
      if (node_bucket_idx < 0 || !m_bucket_activated[node_bucket_idx]) continue;

      const int node_idx = indices_y(i, 1);
      if (node_idx < 0) continue;

      Scalar fnv = m_node_vel_fluid_plus_y[node_bucket_idx](node_idx);

      fv(1) += fnv * weights(i, 1);
    }

    for (int i = 0; i < indices_z.rows(); ++i) {
      const int node_bucket_idx = indices_z(i, 0);
      if (node_bucket_idx < 0 || !m_bucket_activated[node_bucket_idx]) continue;

      const int node_idx = indices_z(i, 1);
      if (node_idx < 0) continue;

      Scalar fnv = m_node_vel_fluid_plus_z[node_bucket_idx](node_idx);

      fv(2) += fnv * weights(i, 2);
    }

    const VecXx& disp = m_strands[strand_idx]->dynamics().getDisplacements();
    const Vec3x u_s = (disp.segment<3>(local_idx * 4) +
                       disp.segment<3>((local_idx + 1) * 4)) *
                      0.5 / m_dt;

    const Vec3x e =
        m_strands[strand_idx]->getEdgeVector(local_idx).normalized();
    const Scalar u_tau = (fv - u_s).dot(e);

    m_strands[strand_idx]->getStepper()->setFutureFlowVelocity(local_idx,
                                                               u_tau);
  });

  // check_isnan("flowvel_after_constraint", (int) m_strands.size(), [this] (int
  // i) {return m_strands[i]->getStepper()->getFutureFlowVelocity();});
}

void LiquidSimulator::updateSurfFlowConstraint() {
  // map flow velocity to grid
  const Scalar dx = getDx();
  const Scalar dV = dx * dx * dx;

  // map interfacing flow to grid
  m_edge_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    const auto& bucket_node_edges_x = m_node_edge_x[bucket_idx];
    const auto& bucket_node_edges_y = m_node_edge_y[bucket_idx];
    const auto& bucket_node_edges_z = m_node_edge_z[bucket_idx];

    VecXx& bucket_node_coeff_x = m_node_vel_constraint_coeff_x[bucket_idx];
    VecXx& bucket_node_coeff_y = m_node_vel_constraint_coeff_y[bucket_idx];
    VecXx& bucket_node_coeff_z = m_node_vel_constraint_coeff_z[bucket_idx];

    bucket_node_coeff_x.setZero();
    bucket_node_coeff_y.setZero();
    bucket_node_coeff_z.setZero();

    VecXx& bucket_node_rhs_x = m_node_rhs_fluid_x[bucket_idx];
    VecXx& bucket_node_rhs_y = m_node_rhs_fluid_y[bucket_idx];
    VecXx& bucket_node_rhs_z = m_node_rhs_fluid_z[bucket_idx];

    const int num_nodes = getNumNodes(bucket_idx);

    for (int i = 0; i < num_nodes; ++i) {
      const auto& node_edge_x = bucket_node_edges_x[i];

      // ignore air nodes
      if (m_node_mass_fluid_x[bucket_idx][i] < 1e-20) continue;

      Scalar p = 0.0;
      Scalar m = 0.0;
      for (auto& pair : node_edge_x) {
        const int pidx = pair.first;
        const int strand_idx = m_strands_inside_segs[pidx].first;
        const int local_idx = m_strands_inside_segs[pidx].second;
        auto& weights = m_inside_edge_weights[pidx];

        const Scalar w = weights(pair.second, 0);

        const Scalar mass =
            m_strands[strand_idx]->dynamics().getFlowMasses()(local_idx);
        const Scalar u_tau =
            m_strands[strand_idx]->getStepper()->getFutureFlowVelocity()(
                local_idx);

        const VecXx& vels = m_strands[strand_idx]->getStepper()->velocities();
        const Scalar u_s =
            (vels(local_idx * 4 + 0) + vels((local_idx + 1) * 4 + 0)) * 0.5;
        const Vec3x e =
            m_strands[strand_idx]->getEdgeVector(local_idx).normalized();

        p += (u_s + e(0) * u_tau) * mass * w;
        m += mass * w;
      }

      bucket_node_rhs_x[i] += p;
      bucket_node_coeff_x[i] = m;
    }

    for (int i = 0; i < num_nodes; ++i) {
      const auto& node_edge_y = bucket_node_edges_y[i];

      // ignore air nodes
      if (m_node_mass_fluid_y[bucket_idx][i] < 1e-20) continue;

      Scalar p = 0.0;
      Scalar m = 0.0;
      for (auto& pair : node_edge_y) {
        const int pidx = pair.first;
        const int strand_idx = m_strands_inside_segs[pidx].first;
        const int local_idx = m_strands_inside_segs[pidx].second;
        auto& weights = m_inside_edge_weights[pidx];

        const Scalar w = weights(pair.second, 1);

        const Scalar mass =
            m_strands[strand_idx]->dynamics().getFlowMasses()(local_idx);
        const Scalar u_tau =
            m_strands[strand_idx]->getStepper()->getFutureFlowVelocity()(
                local_idx);

        const VecXx& vels = m_strands[strand_idx]->getStepper()->velocities();
        const Scalar u_s =
            (vels(local_idx * 4 + 1) + vels((local_idx + 1) * 4 + 1)) * 0.5;
        const Vec3x e =
            m_strands[strand_idx]->getEdgeVector(local_idx).normalized();

        p += (u_s + e(1) * u_tau) * mass * w;
        m += mass * w;
      }

      bucket_node_rhs_y[i] += p;
      bucket_node_coeff_y[i] = m;
    }

    for (int i = 0; i < num_nodes; ++i) {
      const auto& node_edge_z = bucket_node_edges_z[i];

      // ignore air nodes
      if (m_node_mass_fluid_z[bucket_idx][i] < 1e-20) continue;

      Scalar p = 0.0;
      Scalar m = 0.0;
      for (auto& pair : node_edge_z) {
        const int pidx = pair.first;
        const int strand_idx = m_strands_inside_segs[pidx].first;
        const int local_idx = m_strands_inside_segs[pidx].second;
        auto& weights = m_inside_edge_weights[pidx];

        const Scalar w = weights(pair.second, 2);

        const Scalar mass =
            m_strands[strand_idx]->dynamics().getFlowMasses()(local_idx);
        const Scalar u_tau =
            m_strands[strand_idx]->getStepper()->getFutureFlowVelocity()(
                local_idx);

        const VecXx& vels = m_strands[strand_idx]->getStepper()->velocities();
        const Scalar u_s =
            (vels(local_idx * 4 + 2) + vels((local_idx + 1) * 4 + 2)) * 0.5;
        const Vec3x e =
            m_strands[strand_idx]->getEdgeVector(local_idx).normalized();

        p += (u_s + e(2) * u_tau) * mass * w;
        m += mass * w;
      }

      bucket_node_rhs_z[i] += p;
      bucket_node_coeff_z[i] = m;
    }
  });

  // check_isnan("rhs_after_constraint", m_node_rhs_fluid_x, m_node_rhs_fluid_y,
  // m_node_rhs_fluid_z); check_isnan("coeff_after_constraint",
  // m_node_vel_constraint_coeff_x, m_node_vel_constraint_coeff_y,
  // m_node_vel_constraint_coeff_z);
}

void LiquidSimulator::computeRHS() {
  const Scalar dx = getDx();
  const Scalar dV = dx * dx * dx;
  const Scalar iD = getInverseDCoeff();

  const Vec3x g = GravitationForce::getGravity();

  // check_isnan("mass before compute rhs", m_node_mass_fluid_x,
  // m_node_mass_fluid_y, m_node_mass_fluid_z); check_isnan("drag before compute
  // rhs", m_node_drag_coeff_x, m_node_drag_coeff_y, m_node_drag_coeff_z);
  // check_isnan("vel elastic before compute rhs", m_node_vel_elastic_x,
  // m_node_vel_elastic_y, m_node_vel_elastic_z); check_isnan("vel constraint
  // before compute rhs", m_node_vel_constraint_coeff_x,
  // m_node_vel_constraint_coeff_y, m_node_vel_constraint_coeff_z);
  // check_isnan("vol before compute rhs", m_vol);
  // check_isnan("x before compute rhs", m_x);

  if (!m_liquid_info.use_implicit_pressure) {
    m_particle_buckets.for_each_bucket([&](int bucket_idx) {
      if (!m_bucket_activated[bucket_idx]) return;

      const auto& bucket_node_particles_x = m_node_particles_x[bucket_idx];
      const auto& bucket_node_particles_y = m_node_particles_y[bucket_idx];
      const auto& bucket_node_particles_z = m_node_particles_z[bucket_idx];
      const auto& bucket_node_particles_p = m_node_particles_p[bucket_idx];

      const auto& bucket_node_vertex_x = m_node_vertex_x[bucket_idx];
      const auto& bucket_node_vertex_y = m_node_vertex_y[bucket_idx];
      const auto& bucket_node_vertex_z = m_node_vertex_z[bucket_idx];
      const auto& bucket_node_vertex_p = m_node_vertex_p[bucket_idx];

      const int num_nodes = getNumNodes(bucket_idx);

      for (int i = 0; i < num_nodes; ++i) {
        const auto& node_particles_x = bucket_node_particles_x[i];
        const Vec3x& np = getNodePosX(bucket_idx, i);

        Scalar force = m_node_mass_fluid_x[bucket_idx][i] * g(0) +
                       m_node_drag_coeff_x[bucket_idx][i] *
                           m_node_vel_elastic_x[bucket_idx][i];

        if (m_liquid_info.solve_viscosity) {
          for (auto& pair : node_particles_x) {
            const int pidx = pair.first;

            //                    auto& grad_x = m_particle_grad_x[pidx];
            const Scalar v0 = m_vol(pidx);
            const Vec3x& pos = m_x.segment<3>(pidx * 3);
            const Mat3x devb = m_b_plus.block<3, 3>(pidx * 3, 0) -
                               m_b_plus.block<3, 3>(pidx * 3, 0).trace() / 3.0 *
                                   Mat3x::Identity();
            const VecXx& color =
                m_components.segment(pidx * m_liquid_info.num_components,
                                     m_liquid_info.num_components);
            const Scalar pressure =
                getBulkModulus(color) * 0.5 * (m_J(pidx) * m_J(pidx) - 1.0);

            // check_isnan("devb during rhs x", devb.sum());
            // check_isnan("w during rhs x",
            // m_particle_weights[pidx](pair.second, 0));

            const Scalar liquid_shear_modulus = getShearModulus(color);
            const Vec3x dp = np - pos;

            force += -(liquid_shear_modulus * devb.row(0).dot(dp) +
                       pressure * dp(0)) *
                     v0 * iD * m_particle_weights[pidx](pair.second, 0);
          }
        }

        m_node_rhs_fluid_x[bucket_idx](i) += force * m_dt;
      }

      for (int i = 0; i < num_nodes; ++i) {
        const auto& node_particles_y = bucket_node_particles_y[i];
        const Vec3x& np = getNodePosY(bucket_idx, i);

        Scalar force = m_node_mass_fluid_y[bucket_idx][i] * g(1) +
                       m_node_drag_coeff_y[bucket_idx][i] *
                           m_node_vel_elastic_y[bucket_idx][i];

        if (m_liquid_info.solve_viscosity) {
          for (auto& pair : node_particles_y) {
            const int pidx = pair.first;

            //                    auto& grad_y = m_particle_grad_y[pidx];
            const Scalar v0 = m_vol(pidx);
            const Vec3x& pos = m_x.segment<3>(pidx * 3);
            const Mat3x devb = m_b_plus.block<3, 3>(pidx * 3, 0) -
                               m_b_plus.block<3, 3>(pidx * 3, 0).trace() / 3.0 *
                                   Mat3x::Identity();
            const VecXx& color =
                m_components.segment(pidx * m_liquid_info.num_components,
                                     m_liquid_info.num_components);
            const Scalar pressure =
                getBulkModulus(color) * 0.5 * (m_J(pidx) * m_J(pidx) - 1.0);

            // check_isnan("devb during rhs y", devb.sum());
            // check_isnan("w during rhs y",
            // m_particle_weights[pidx](pair.second, 1));

            const Scalar liquid_shear_modulus = getShearModulus(color);
            const Vec3x dp = np - pos;

            force += -(liquid_shear_modulus * devb.row(1).dot(dp) +
                       pressure * dp(1)) *
                     v0 * iD * m_particle_weights[pidx](pair.second, 1);
          }
        }

        m_node_rhs_fluid_y[bucket_idx](i) += force * m_dt;
      }

      for (int i = 0; i < num_nodes; ++i) {
        const auto& node_particles_z = bucket_node_particles_z[i];
        const Vec3x& np = getNodePosZ(bucket_idx, i);

        Scalar force = m_node_mass_fluid_z[bucket_idx][i] * g(2) +
                       m_node_drag_coeff_z[bucket_idx][i] *
                           m_node_vel_elastic_z[bucket_idx][i];

        if (m_liquid_info.solve_viscosity) {
          for (auto& pair : node_particles_z) {
            const int pidx = pair.first;

            //                    auto& grad_z = m_particle_grad_z[pidx];
            const Scalar v0 = m_vol(pidx);
            const Vec3x& pos = m_x.segment<3>(pidx * 3);
            const Mat3x devb = m_b_plus.block<3, 3>(pidx * 3, 0) -
                               m_b_plus.block<3, 3>(pidx * 3, 0).trace() / 3.0 *
                                   Mat3x::Identity();
            const VecXx& color =
                m_components.segment(pidx * m_liquid_info.num_components,
                                     m_liquid_info.num_components);
            const Scalar pressure =
                getBulkModulus(color) * 0.5 * (m_J(pidx) * m_J(pidx) - 1.0);

            // check_isnan("devb during rhs z", devb.sum());
            // check_isnan("w during rhs z",
            // m_particle_weights[pidx](pair.second, 2));

            const Scalar liquid_shear_modulus = getShearModulus(color);
            const Vec3x dp = np - pos;

            force += -(liquid_shear_modulus * devb.row(2).dot(dp) +
                       pressure * dp(2)) *
                     v0 * iD * m_particle_weights[pidx](pair.second, 2);
          }
        }

        m_node_rhs_fluid_z[bucket_idx](i) += force * m_dt;
      }
    });
  } else {
    m_particle_buckets.for_each_bucket([&](int bucket_idx) {
      if (!m_bucket_activated[bucket_idx]) return;

      const auto& bucket_node_particles_x = m_node_particles_x[bucket_idx];
      const auto& bucket_node_particles_y = m_node_particles_y[bucket_idx];
      const auto& bucket_node_particles_z = m_node_particles_z[bucket_idx];
      const auto& bucket_node_particles_p = m_node_particles_p[bucket_idx];

      const auto& bucket_node_vertex_x = m_node_vertex_x[bucket_idx];
      const auto& bucket_node_vertex_y = m_node_vertex_y[bucket_idx];
      const auto& bucket_node_vertex_z = m_node_vertex_z[bucket_idx];
      const auto& bucket_node_vertex_p = m_node_vertex_p[bucket_idx];

      const int num_nodes = getNumNodes(bucket_idx);

      for (int i = 0; i < num_nodes; ++i) {
        const auto& node_particles_x = bucket_node_particles_x[i];
        const Vec3x& np = getNodePosX(bucket_idx, i);

        Scalar force = m_node_mass_fluid_x[bucket_idx][i] * g(0) +
                       m_node_drag_coeff_x[bucket_idx][i] *
                           m_node_vel_elastic_x[bucket_idx][i];

        if (m_liquid_info.solve_viscosity) {
          for (auto& pair : node_particles_x) {
            const int pidx = pair.first;

            const VecXx& color =
                m_components.segment(pidx * m_liquid_info.num_components,
                                     m_liquid_info.num_components);
            const Scalar liquid_shear_modulus = getShearModulus(color);
            if (liquid_shear_modulus > 1e-20) {
              //                    auto& grad_x = m_particle_grad_x[pidx];
              const Scalar v0 = m_vol(pidx);
              const Vec3x& pos = m_x.segment<3>(pidx * 3);
              const Mat3x devb = m_b_plus.block<3, 3>(pidx * 3, 0) -
                                 m_b_plus.block<3, 3>(pidx * 3, 0).trace() /
                                     3.0 * Mat3x::Identity();

              // check_isnan("devb during rhs x", devb.sum());
              // check_isnan("w during rhs x",
              // m_particle_weights[pidx](pair.second, 0));

              force += -liquid_shear_modulus * v0 * devb.row(0).dot(np - pos) *
                       iD * m_particle_weights[pidx](pair.second, 0);
            }
          }
        }

        m_node_rhs_fluid_x[bucket_idx](i) += force * m_dt;
      }

      for (int i = 0; i < num_nodes; ++i) {
        const auto& node_particles_y = bucket_node_particles_y[i];
        const Vec3x& np = getNodePosY(bucket_idx, i);

        Scalar force = m_node_mass_fluid_y[bucket_idx][i] * g(1) +
                       m_node_drag_coeff_y[bucket_idx][i] *
                           m_node_vel_elastic_y[bucket_idx][i];

        if (m_liquid_info.solve_viscosity) {
          for (auto& pair : node_particles_y) {
            const int pidx = pair.first;
            const VecXx& color =
                m_components.segment(pidx * m_liquid_info.num_components,
                                     m_liquid_info.num_components);
            const Scalar liquid_shear_modulus = getShearModulus(color);

            if (liquid_shear_modulus > 1e-20) {
              //                    auto& grad_y = m_particle_grad_y[pidx];
              const Scalar v0 = m_vol(pidx);
              const Vec3x& pos = m_x.segment<3>(pidx * 3);
              const Mat3x devb = m_b_plus.block<3, 3>(pidx * 3, 0) -
                                 m_b_plus.block<3, 3>(pidx * 3, 0).trace() /
                                     3.0 * Mat3x::Identity();

              // check_isnan("devb during rhs y", devb.sum());
              // check_isnan("w during rhs y",
              // m_particle_weights[pidx](pair.second, 1));

              force += -liquid_shear_modulus * v0 * devb.row(1).dot(np - pos) *
                       iD * m_particle_weights[pidx](pair.second, 1);
            }
          }
        }

        m_node_rhs_fluid_y[bucket_idx](i) += force * m_dt;
      }

      for (int i = 0; i < num_nodes; ++i) {
        const auto& node_particles_z = bucket_node_particles_z[i];
        const Vec3x& np = getNodePosZ(bucket_idx, i);

        Scalar force = m_node_mass_fluid_z[bucket_idx][i] * g(2) +
                       m_node_drag_coeff_z[bucket_idx][i] *
                           m_node_vel_elastic_z[bucket_idx][i];

        if (m_liquid_info.solve_viscosity) {
          for (auto& pair : node_particles_z) {
            const int pidx = pair.first;
            const VecXx& color =
                m_components.segment(pidx * m_liquid_info.num_components,
                                     m_liquid_info.num_components);
            const Scalar liquid_shear_modulus = getShearModulus(color);

            if (liquid_shear_modulus > 1e-20) {
              //                    auto& grad_z = m_particle_grad_z[pidx];
              const Scalar v0 = m_vol(pidx);
              const Vec3x& pos = m_x.segment<3>(pidx * 3);
              const Mat3x devb = m_b_plus.block<3, 3>(pidx * 3, 0) -
                                 m_b_plus.block<3, 3>(pidx * 3, 0).trace() /
                                     3.0 * Mat3x::Identity();

              // check_isnan("devb during rhs z", devb.sum());
              // check_isnan("w during rhs z",
              // m_particle_weights[pidx](pair.second, 2));

              force += -liquid_shear_modulus * v0 * devb.row(2).dot(np - pos) *
                       iD * m_particle_weights[pidx](pair.second, 2);
            }
          }
        }

        m_node_rhs_fluid_z[bucket_idx](i) += force * m_dt;
      }
    });
  }
  // check_isnan("mass after compute rhs", m_node_mass_fluid_x,
  // m_node_mass_fluid_y, m_node_mass_fluid_z); check_isnan("drag after compute
  // rhs", m_node_drag_coeff_x, m_node_drag_coeff_y, m_node_drag_coeff_z);
  // check_isnan("vel elastic after compute rhs", m_node_vel_elastic_x,
  // m_node_vel_elastic_y, m_node_vel_elastic_z); check_isnan("vel constraint
  // after compute rhs", m_node_vel_constraint_coeff_x,
  // m_node_vel_constraint_coeff_y, m_node_vel_constraint_coeff_z);
  // check_isnan("vol after compute rhs", m_vol);
  // check_isnan("x after compute rhs", m_x);
}

Scalar LiquidSimulator::lengthNodeVectors(
    const std::vector<VecXx>& node_vec_ax,
    const std::vector<VecXx>& node_vec_ay,
    const std::vector<VecXx>& node_vec_az) const {
  VecXx bucket_length(node_vec_ax.size());

  const int num_buckets = node_vec_ax.size();

  for_each(0, num_buckets, [&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) {
      bucket_length(bucket_idx) = 0.0;
    } else {
      bucket_length[bucket_idx] = node_vec_ax[bucket_idx].squaredNorm() +
                                  node_vec_ay[bucket_idx].squaredNorm() +
                                  node_vec_az[bucket_idx].squaredNorm();
    }
  });

  return sqrt(bucket_length.sum());
}

Scalar LiquidSimulator::dotNodeVectors(
    const std::vector<VecXx>& node_vec_ax,
    const std::vector<VecXx>& node_vec_ay,
    const std::vector<VecXx>& node_vec_az,
    const std::vector<VecXx>& node_vec_bx,
    const std::vector<VecXx>& node_vec_by,
    const std::vector<VecXx>& node_vec_bz) const {
  VecXx bucket_dot(node_vec_ax.size());

  const int num_buckets = node_vec_ax.size();

  for_each(0, num_buckets, [&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) {
      bucket_dot(bucket_idx) = 0.0;
    } else {
      bucket_dot[bucket_idx] =
          node_vec_ax[bucket_idx].dot(node_vec_bx[bucket_idx]) +
          node_vec_ay[bucket_idx].dot(node_vec_by[bucket_idx]) +
          node_vec_az[bucket_idx].dot(node_vec_bz[bucket_idx]);
    }
  });

  return bucket_dot.sum();
}

Vec3x LiquidSimulator::addForce(const ElasticStrand*, int) {
  return Vec3x::Zero();
}
Mat3x LiquidSimulator::addJacobian(const ElasticStrand*, int) {
  return Mat3x::Zero();
}
bool LiquidSimulator::advectParticles() {
  LockGuard lock(m_particle_mutex);

  m_x += m_v * m_dt;
  return true;
}
void LiquidSimulator::solveForceExplicit(
    const vector<VecXx>& lhs_x, const vector<VecXx>& lhs_y,
    const vector<VecXx>& lhs_z, const vector<VecXx>& rhs_x,
    const vector<VecXx>& rhs_y, const vector<VecXx>& rhs_z,
    vector<VecXx>& ret_x, vector<VecXx>& ret_y, vector<VecXx>& ret_z) {
  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    const int num_nodes = getNumNodes(bucket_idx);

    for (int i = 0; i < num_nodes; ++i) {
      if (lhs_x[bucket_idx][i] > 1e-20) {
        ret_x[bucket_idx][i] = rhs_x[bucket_idx][i] / lhs_x[bucket_idx][i];
      } else {
        ret_x[bucket_idx][i] = 0.0;
      }

      if (lhs_y[bucket_idx][i] > 1e-20) {
        ret_y[bucket_idx][i] = rhs_y[bucket_idx][i] / lhs_y[bucket_idx][i];
      } else {
        ret_y[bucket_idx][i] = 0.0;
      }

      if (lhs_z[bucket_idx][i] > 1e-20) {
        ret_z[bucket_idx][i] = rhs_z[bucket_idx][i] / lhs_z[bucket_idx][i];
      } else {
        ret_z[bucket_idx][i] = 0.0;
      }
    }
  });
}
void LiquidSimulator::solveForceExplicit(const vector<VecXx>& lhs_x,
                                         const vector<VecXx>& lhs_y,
                                         const vector<VecXx>& lhs_z) {
  solveForceExplicit(lhs_x, lhs_y, lhs_z, m_node_rhs_fluid_x,
                     m_node_rhs_fluid_y, m_node_rhs_fluid_z,
                     m_node_vel_fluid_plus_x, m_node_vel_fluid_plus_y,
                     m_node_vel_fluid_plus_z);
}

bool LiquidSimulator::addForceFluids() {
  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    // init next step velocity
    m_node_vel_fluid_plus_x[bucket_idx] = m_node_vel_fluid_x[bucket_idx];
    m_node_vel_fluid_plus_y[bucket_idx] = m_node_vel_fluid_y[bucket_idx];
    m_node_vel_fluid_plus_z[bucket_idx] = m_node_vel_fluid_z[bucket_idx];

    m_node_rhs_fluid_x[bucket_idx] =
        VecXx(m_node_mass_fluid_x[bucket_idx].array() *
              m_node_vel_fluid_x[bucket_idx].array());
    m_node_rhs_fluid_y[bucket_idx] =
        VecXx(m_node_mass_fluid_y[bucket_idx].array() *
              m_node_vel_fluid_y[bucket_idx].array());
    m_node_rhs_fluid_z[bucket_idx] =
        VecXx(m_node_mass_fluid_z[bucket_idx].array() *
              m_node_vel_fluid_z[bucket_idx].array());
  });

  applyGridBasedBodyFriction(m_node_mass_fluid_x, m_node_mass_fluid_y,
                             m_node_mass_fluid_z);

  // check_isnan("fluidvel_after_friction", m_node_vel_fluid_plus_x,
  // m_node_vel_fluid_plus_y, m_node_vel_fluid_plus_z);

  findInterfacingSegments();
  updateSurfFlowConstraint();

  computeRHS();

  // check_isnan("rhs_after_computeRHS", m_node_rhs_fluid_x, m_node_rhs_fluid_y,
  // m_node_rhs_fluid_z);

  computeShearLHSDiagonal();

  // check_isnan("lhs_after_computeLHS", m_node_lhs_fluid_x, m_node_lhs_fluid_y,
  // m_node_lhs_fluid_z);

  return true;
}

Vec3x LiquidSimulator::get_previous_velocity(const Vec3x& position) const {
  return Vec3x::Zero();
}

bool LiquidSimulator::saveVelocity() { return true; }

void LiquidSimulator::solveColorDiffusion() {
  if (!m_liquid_info.solve_color_diffusion) return;

  const int num_fluid = numParticles();
  const Scalar dx = getDx();

  m_particle_cells_single.sort(
      num_fluid, [&](int pidx, int& i, int& j, int& k) {
        Vec3x local_x = (m_x.segment<3>(pidx * 3) - m_grid_mincorner) / dx;
        i = (int)floor(local_x(0));
        j = (int)floor(local_x(1));
        k = (int)floor(local_x(2));
      });

  VecXx sph_density(num_fluid);
  VecXx mu(num_fluid * m_liquid_info.num_components);

  sph_density.setZero();
  mu.setZero();

  const Scalar inv_dx3 = 1.0 / (dx * dx * dx);

  m_particle_cells_single.for_each_bucket_particles(
      [&](int pidx, int bucket_idx) {
        Scalar sum_w = 0.0;
        const Vec3x& pos = m_x.segment<3>(pidx * 3);

        m_particle_cells_single.loop_neighbor_bucket_particles(
            bucket_idx, [&](int npidx, int) {
              const Vec3x& npos = m_x.segment<3>(npidx * 3);

              const Scalar r2 = (pos - npos).squaredNorm();
              const Scalar w = smooth_kernel(r2, dx);
              sum_w += w;

              return false;
            });

        sph_density(pidx) = m_m(pidx * 3) * sum_w * 1.5666814711 * inv_dx3;
      });

  // check_isnan("sph_density", sph_density);

  const Scalar eta2 = 0.01 * dx * dx;
  const Scalar inv_num_components = 1.0 / (Scalar)m_liquid_info.num_components;
  const Scalar epsilon2 =
      m_liquid_info.geometric_diffusivity * m_liquid_info.geometric_diffusivity;

  m_particle_cells_single.for_each_bucket_particles(
      [&](int pidx, int bucket_idx) {
        VecXx local_mu = VecXx::Zero(m_liquid_info.num_components);

        const Vec3x& pos = m_x.segment<3>(pidx * 3);
        const VecXx& color = m_components.segment(
            pidx * m_liquid_info.num_components, m_liquid_info.num_components);

        m_particle_cells_single.loop_neighbor_bucket_particles(
            bucket_idx, [&](int npidx, int) {
              const Vec3x& npos = m_x.segment<3>(npidx * 3);

              const Scalar r2 = (pos - npos).squaredNorm();

              const Scalar gwr = grad_smooth_kernel_dot_dist(r2, dx);

              const Scalar coeff = m_m(npidx * 3) /
                                   std::max(1e-20, sph_density(npidx)) * gwr /
                                   (r2 + eta2);

              local_mu += (color - m_components.segment(
                                       npidx * m_liquid_info.num_components,
                                       m_liquid_info.num_components)) *
                          coeff;

              return false;
            });

        local_mu *= 18.8001776528 * inv_dx3 * epsilon2;

        Scalar beta_c = 0.0;

        for (int k = 0; k < m_liquid_info.num_components; ++k) {
          Scalar dfdck = m_liquid_info.chemical_diffusivity *
                         (color(k) - m_liquid_info.helmholtz_first_minima(k)) *
                         2. * inv_num_components;

          beta_c -= dfdck;
          local_mu(k) += dfdck;
        }

        beta_c *= inv_num_components;

        local_mu += VecXx::Constant(m_liquid_info.num_components, beta_c);

        mu.segment(pidx * m_liquid_info.num_components,
                   m_liquid_info.num_components) = local_mu;
      });

  // check_isnan("mu", mu);

  m_particle_cells_single.for_each_bucket_particles(
      [&](int pidx, int bucket_idx) {
        VecXx local_laplacian_mu = VecXx::Zero(m_liquid_info.num_components);

        const Vec3x& pos = m_x.segment<3>(pidx * 3);

        const VecXx& mu_i = mu.segment(pidx * m_liquid_info.num_components,
                                       m_liquid_info.num_components);

        m_particle_cells_single.loop_neighbor_bucket_particles(
            bucket_idx, [&](int npidx, int) {
              const Vec3x& npos = m_x.segment<3>(npidx * 3);

              const Scalar r2 = (pos - npos).squaredNorm();

              const Scalar gwr = grad_smooth_kernel_dot_dist(r2, dx);

              const Scalar coeff = m_m(npidx * 3) /
                                   std::max(1e-20, sph_density(npidx)) * gwr /
                                   (r2 + eta2);

              local_laplacian_mu +=
                  (mu_i - mu.segment(npidx * m_liquid_info.num_components,
                                     m_liquid_info.num_components)) *
                  coeff;

              return false;
            });

        local_laplacian_mu *= -18.8001776528 * inv_dx3;

        VecXx color = m_components.segment(pidx * m_liquid_info.num_components,
                                           m_liquid_info.num_components) +
                      local_laplacian_mu * m_dt;

        for (int k = 0; k < m_liquid_info.num_components; ++k) {
          color(k) = clamp(color(k), 0.0, 1.0);
        }

        make_gibbs_simplex(color);

        m_components.segment(pidx * m_liquid_info.num_components,
                             m_liquid_info.num_components) = color;
      });

  // check_isnan("components", m_components);

  //        std::cout << m_components << std::endl;
}

void LiquidSimulator::multiplyHessianMatrix(const std::vector<VecXx>& vec_x,
                                            const std::vector<VecXx>& vec_y,
                                            const std::vector<VecXx>& vec_z,
                                            std::vector<VecXx>& ret_x,
                                            std::vector<VecXx>& ret_y,
                                            std::vector<VecXx>& ret_z) {
  if (!m_liquid_info.solve_viscosity) {
    m_particle_buckets.for_each_bucket([&](int bucket_idx) {
      if (!m_bucket_activated[bucket_idx]) return;

      const VecXx& bucket_vec_x = vec_x[bucket_idx];
      const VecXx& bucket_vec_y = vec_y[bucket_idx];
      const VecXx& bucket_vec_z = vec_z[bucket_idx];

      ret_x[bucket_idx] = VecXx((m_node_mass_fluid_x[bucket_idx] +
                                 m_node_drag_coeff_x[bucket_idx] * m_dt +
                                 m_node_vel_constraint_coeff_x[bucket_idx])
                                    .array() *
                                bucket_vec_x.array());
      ret_y[bucket_idx] = VecXx((m_node_mass_fluid_y[bucket_idx] +
                                 m_node_drag_coeff_y[bucket_idx] * m_dt +
                                 m_node_vel_constraint_coeff_y[bucket_idx])
                                    .array() *
                                bucket_vec_y.array());
      ret_z[bucket_idx] = VecXx((m_node_mass_fluid_z[bucket_idx] +
                                 m_node_drag_coeff_z[bucket_idx] * m_dt +
                                 m_node_vel_constraint_coeff_z[bucket_idx])
                                    .array() *
                                bucket_vec_z.array());
    });
  } else {
    const int num_parts = numParticles();
    m_Ap.resize(num_parts * 3, 3);

    const Scalar iD = getInverseDCoeff();

    auto flat_mul = [](const Mat3x& a, const Mat3x& b) -> Scalar {
      return a.col(0).dot(b.col(0)) + a.col(1).dot(b.col(1)) +
             a.col(2).dot(b.col(2));
    };

    // Calculate A_p and DW for each particle
    for_each(0, num_parts, [&](int pidx) {
      const VecXx& color = m_components.segment(
          m_liquid_info.num_components * pidx, m_liquid_info.num_components);

      const Scalar mu = getShearModulus(color);

      if (mu > 1e-20) {
        auto& indices_x = m_particle_nodes_x[pidx];
        auto& indices_y = m_particle_nodes_y[pidx];
        auto& indices_z = m_particle_nodes_z[pidx];

        const Vec3x& pos = m_x.segment<3>(pidx * 3);

        Mat3x B = Mat3x::Zero();

        const Mat27x4f& weights = m_particle_weights[pidx];

        auto add_to_B = [&](const vector<VecXx>& vec, const Mat27x2i& indices,
                            int ir) {
          for (int i = 0; i < indices.rows(); ++i) {
            const int node_bucket_idx = indices(i, 0);
            if (!m_bucket_activated[node_bucket_idx]) {
              continue;
            }

            const int node_idx = indices(i, 1);
            Vec3x np = getNodePos(node_bucket_idx, node_idx, ir);

            B.col(ir) += (np - pos) * iD * weights(i, ir) *
                         vec[node_bucket_idx][node_idx];
          }
        };

        add_to_B(vec_x, indices_x, 0);
        add_to_B(vec_y, indices_y, 1);
        add_to_B(vec_z, indices_z, 2);

        const Mat3x& b_bar = m_b_trial.block<3, 3>(pidx * 3, 0);
        const Scalar Ie_bar = b_bar.trace() / 3.0;

        const Mat3x devb_bar = b_bar - Ie_bar * Mat3x::Identity();

        const Scalar devb_bar_n = devb_bar.norm();
        const Scalar BTb_bar_tr = (B.transpose() * b_bar).trace();
        const Scalar b_bar_tr = b_bar.trace();

        // Elastic Part
        Mat3x Ap =
            (B.transpose() * b_bar - 2.0 / 3.0 * B.trace() * devb_bar -
             1.0 / 3.0 * (BTb_bar_tr * Mat3x::Identity() - b_bar_tr * B)) *
            m_proj_func(pidx * 3);

        const Scalar s = mu * devb_bar_n;
        const Scalar tilde_sigma_Y = sqrt(2.0 / 3.0) * getYieldStress(color);
        const Scalar phi_trial = s - tilde_sigma_Y;

        if (s > 0.0) {
          const Scalar dgds = m_proj_func(pidx * 3 + 1);
          const Scalar dgdtmu = m_proj_func(pidx * 3 + 2);
          const Scalar B_tr = B.trace();

          // Plastic Part
          Mat3x plastic_part =
              devb_bar *
              (dgds * (mu * (devb_bar * B.transpose() * b_bar).trace() /
                           devb_bar_n -
                       2.0 * s / 3.0 * B_tr) +
               dgdtmu * mu / 3.0 * (BTb_bar_tr - 2.0 / 3.0 * B_tr * b_bar_tr));

          Ap += plastic_part;
        }

        const Scalar coeff = mu * m_dt * m_dt;

        m_Ap.block<3, 3>(pidx * 3, 0) =
            m_vol(pidx) * coeff * Ap *
            (1.0 + m_liquid_info.liquid_shear_damping / m_dt);
      } else {
        m_Ap.block<3, 3>(pidx * 3, 0).setZero();
      }
    });

    // Weighted Sum of A_p to results
    m_particle_buckets.for_each_bucket([&](int bucket_idx) {
      if (!m_bucket_activated[bucket_idx]) return;

      const VecXx& bucket_vec_x = vec_x[bucket_idx];
      const VecXx& bucket_vec_y = vec_y[bucket_idx];
      const VecXx& bucket_vec_z = vec_z[bucket_idx];

      ret_x[bucket_idx] = VecXx((m_node_mass_fluid_x[bucket_idx] +
                                 m_node_drag_coeff_x[bucket_idx] * m_dt +
                                 m_node_vel_constraint_coeff_x[bucket_idx])
                                    .array() *
                                bucket_vec_x.array());
      ret_y[bucket_idx] = VecXx((m_node_mass_fluid_y[bucket_idx] +
                                 m_node_drag_coeff_y[bucket_idx] * m_dt +
                                 m_node_vel_constraint_coeff_y[bucket_idx])
                                    .array() *
                                bucket_vec_y.array());
      ret_z[bucket_idx] = VecXx((m_node_mass_fluid_z[bucket_idx] +
                                 m_node_drag_coeff_z[bucket_idx] * m_dt +
                                 m_node_vel_constraint_coeff_z[bucket_idx])
                                    .array() *
                                bucket_vec_z.array());

      const auto& bucket_node_particles_x = m_node_particles_x[bucket_idx];
      const auto& bucket_node_particles_y = m_node_particles_y[bucket_idx];
      const auto& bucket_node_particles_z = m_node_particles_z[bucket_idx];

      const auto& bucket_node_vertex_x = m_node_vertex_x[bucket_idx];
      const auto& bucket_node_vertex_y = m_node_vertex_y[bucket_idx];
      const auto& bucket_node_vertex_z = m_node_vertex_z[bucket_idx];

      const int num_nodes = getNumNodes(bucket_idx);

      auto add_Ap_to_ret = [&](const vector<pair<int, short> >& node_particles,
                               int i, int ir) -> Scalar {
        Scalar ret = 0.0;
        const Vec3x& np = getNodePos(bucket_idx, i, ir);

        for (auto& pair : node_particles) {
          const int pidx = pair.first;
          const Vec3x& pos = m_x.segment<3>(pidx * 3);

          ret += m_Ap.block<1, 3>(pidx * 3 + ir, 0).dot(np - pos) *
                 m_particle_weights[pidx](pair.second, ir) * iD;
        }
        return ret;
      };

      for (int i = 0; i < num_nodes; ++i) {
        ret_x[bucket_idx][i] += add_Ap_to_ret(bucket_node_particles_x[i], i, 0);
        ret_y[bucket_idx][i] += add_Ap_to_ret(bucket_node_particles_y[i], i, 1);
        ret_z[bucket_idx][i] += add_Ap_to_ret(bucket_node_particles_z[i], i, 2);
      }
    });
  }
}

void LiquidSimulator::computeShearLHSDiagonal() {
  allocateNodeVectors(m_node_lhs_fluid_x);
  allocateNodeVectors(m_node_lhs_fluid_y);
  allocateNodeVectors(m_node_lhs_fluid_z);

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    m_node_lhs_fluid_x[bucket_idx] = m_node_mass_fluid_x[bucket_idx] +
                                     m_node_drag_coeff_x[bucket_idx] * m_dt +
                                     m_node_vel_constraint_coeff_x[bucket_idx];
    m_node_lhs_fluid_y[bucket_idx] = m_node_mass_fluid_y[bucket_idx] +
                                     m_node_drag_coeff_y[bucket_idx] * m_dt +
                                     m_node_vel_constraint_coeff_y[bucket_idx];
    m_node_lhs_fluid_z[bucket_idx] = m_node_mass_fluid_z[bucket_idx] +
                                     m_node_drag_coeff_z[bucket_idx] * m_dt +
                                     m_node_vel_constraint_coeff_z[bucket_idx];
  });

  // check_isnan("lhs_diagonal_mass", m_node_mass_fluid_x, m_node_mass_fluid_y,
  // m_node_mass_fluid_z); check_isnan("lhs_diagonal_vel_constraint",
  // m_node_vel_constraint_coeff_x, m_node_vel_constraint_coeff_y,
  // m_node_vel_constraint_coeff_z); check_isnan("lhs_diagonal_drag",
  // m_node_drag_coeff_x, m_node_drag_coeff_y, m_node_drag_coeff_z);
  // check_isnan("lhs_diagonal_lhs", m_node_lhs_fluid_x, m_node_lhs_fluid_y,
  // m_node_lhs_fluid_z);
}

bool LiquidSimulator::computeWeight() {
  //        std::cout << "[compute weight 0]" << std::endl;

  updateParticleWeights();

  //        std::cout << "[compute weight 1]" << std::endl;

  buildNodeParticlePairs();

  //       std::cout << "[compute weight 2]" << std::endl;

  updateVertexWeights();

  //        std::cout << "[compute weight 3]" << std::endl;

  buildNodeVertexPairs();

  //        std::cout << "[compute weight 4]" << std::endl;

  return true;
}
bool LiquidSimulator::computePhi() {
  updateLiquidPhi();
  updateSolidWeights();

  return true;
}

template <typename T>
void LiquidSimulator::allocateNodeVectors(
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& node_vec_x,
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& node_vec_y,
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& node_vec_z) const {
  const int num_buckets = m_particle_buckets.size();

  if ((int)node_vec_x.size() != num_buckets) node_vec_x.resize(num_buckets);
  if ((int)node_vec_y.size() != num_buckets) node_vec_y.resize(num_buckets);
  if ((int)node_vec_z.size() != num_buckets) node_vec_z.resize(num_buckets);

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    const int num_nodes = getNumNodes(bucket_idx);

    node_vec_x[bucket_idx].resize(num_nodes);
    node_vec_x[bucket_idx].setZero();
    node_vec_y[bucket_idx].resize(num_nodes);
    node_vec_y[bucket_idx].setZero();
    node_vec_z[bucket_idx].resize(num_nodes);
    node_vec_z[bucket_idx].setZero();
  });
}

bool LiquidSimulator::doUseDrag() { return m_liquid_info.use_liquid_drag; }

template <typename T>
void LiquidSimulator::allocateNodeVectors(
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& node_vec_p, int num_subs,
    const T& val) const {
  const int num_buckets = m_particle_buckets.size();

  if ((int)node_vec_p.size() != num_buckets) node_vec_p.resize(num_buckets);

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    const int num_nodes = getNumNodes(bucket_idx);

    node_vec_p[bucket_idx].resize(num_nodes * num_subs);

    if (num_nodes > 0) node_vec_p[bucket_idx].setConstant(val);
  });
}

void LiquidSimulator::computeElasticLiquidDifference() {
  const int num_verts = m_strands_global_local.size();

  for_each(0, num_verts, [&](int pidx) {
    const int strand_idx = m_strands_global_local[pidx].first;
    const int local_idx = m_strands_global_local[pidx].second;

    Vec3x elastic_vel =
        m_strands[strand_idx]->dynamics().getDisplacement(local_idx) / m_dt;

    auto& indices_x = m_vertex_nodes_x[pidx];
    auto& indices_y = m_vertex_nodes_y[pidx];
    auto& indices_z = m_vertex_nodes_z[pidx];

    auto& weights = m_vertex_weights[pidx];

    Vec3x fv = Vec3x::Zero();

    for (int i = 0; i < indices_x.rows(); ++i) {
      const int node_bucket_idx = indices_x(i, 0);
      if (node_bucket_idx == -1 || !m_bucket_activated[node_bucket_idx])
        continue;
      const int node_idx = indices_x(i, 1);
      if (node_idx == -1) continue;
      Scalar fnv = m_node_vel_fluid_plus_x[node_bucket_idx](node_idx);
      fv(0) += fnv * weights(i, 0);
    }

    for (int i = 0; i < indices_y.rows(); ++i) {
      const int node_bucket_idx = indices_y(i, 0);
      if (node_bucket_idx == -1 || !m_bucket_activated[node_bucket_idx])
        continue;
      const int node_idx = indices_y(i, 1);
      if (node_idx == -1) continue;
      Scalar fnv = m_node_vel_fluid_plus_y[node_bucket_idx](node_idx);
      fv(1) += fnv * weights(i, 1);
    }

    for (int i = 0; i < indices_z.rows(); ++i) {
      const int node_bucket_idx = indices_z(i, 0);
      if (node_bucket_idx == -1 || !m_bucket_activated[node_bucket_idx])
        continue;
      const int node_idx = indices_z(i, 1);
      if (node_idx == -1) continue;
      Scalar fnv = m_node_vel_fluid_plus_z[node_bucket_idx](node_idx);
      fv(2) += fnv * weights(i, 2);
    }

    const Scalar len = m_strands[strand_idx]->getVoronoiLength(local_idx);
    const Scalar radius = m_strands[strand_idx]->getRadiusA(local_idx);
    const Scalar vol = M_PI * radius * radius * len;

    // check_isnan("eld_fvx", fv(0));
    // check_isnan("eld_fvy", fv(1));
    // check_isnan("eld_fvz", fv(2));

    m_elastic_liquid_diffs.segment<3>(pidx * 3) = vol * (elastic_vel - fv);
  });

  // check_isnan("elastic_liquid_diffs", m_elastic_liquid_diffs);
}

bool LiquidSimulator::solvePressure() {
  if (!m_liquid_info.use_implicit_pressure) return true;

  std::cout << "[compute elastic liquid diff]" << std::endl;
  computeElasticLiquidDifference();

  //        std::cout << "[solvePressure: 0]" << std::endl;
  const int ni = m_particle_buckets.ni * m_num_nodes;
  const int nj = m_particle_buckets.nj * m_num_nodes;
  const int nk = m_particle_buckets.nk * m_num_nodes;

  allocateNodeVectors(m_node_global_pressure_indices);

  std::vector<int> num_effective_nodes(m_particle_buckets.size());
  // assign global indices to nodes
  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    const int num_nodes_p = m_node_liquid_phi[bucket_idx].size();
    const VecXi& bucket_pn = m_node_pressure_neighbors[bucket_idx];

    int k = 0;
    for (int i = 0; i < num_nodes_p; ++i) {
      m_node_global_pressure_indices[bucket_idx][i] = -1;

      if (m_node_liquid_phi[bucket_idx][i] < 0) {
        const int bucket_idx_left = bucket_pn[i * 12 + 0];
        const int node_idx_left = bucket_pn[i * 12 + 1];
        const int bucket_idx_right = bucket_pn[i * 12 + 2];
        const int node_idx_right = bucket_pn[i * 12 + 3];
        const int bucket_idx_bottom = bucket_pn[i * 12 + 4];
        const int node_idx_bottom = bucket_pn[i * 12 + 5];
        const int bucket_idx_top = bucket_pn[i * 12 + 6];
        const int node_idx_top = bucket_pn[i * 12 + 7];
        const int bucket_idx_near = bucket_pn[i * 12 + 8];
        const int node_idx_near = bucket_pn[i * 12 + 9];
        const int bucket_idx_far = bucket_pn[i * 12 + 10];
        const int node_idx_far = bucket_pn[i * 12 + 11];

        if ((bucket_idx_left >= 0 && node_idx_left >= 0 &&
             m_bucket_activated[bucket_idx_left] &&
             m_node_liquid_weight_x[bucket_idx_left][node_idx_left]) ||
            (bucket_idx_right >= 0 && node_idx_right >= 0 &&
             m_bucket_activated[bucket_idx_right] &&
             m_node_liquid_weight_x[bucket_idx_right][node_idx_right]) ||
            (bucket_idx_bottom >= 0 && node_idx_bottom >= 0 &&
             m_bucket_activated[bucket_idx_bottom] &&
             m_node_liquid_weight_y[bucket_idx_bottom][node_idx_bottom]) ||
            (bucket_idx_top >= 0 && node_idx_top >= 0 &&
             m_bucket_activated[bucket_idx_top] &&
             m_node_liquid_weight_y[bucket_idx_top][node_idx_top]) ||
            (bucket_idx_near >= 0 && node_idx_near >= 0 &&
             m_bucket_activated[bucket_idx_near] &&
             m_node_liquid_weight_z[bucket_idx_near][node_idx_near]) ||
            (bucket_idx_far >= 0 && node_idx_far >= 0 &&
             m_bucket_activated[bucket_idx_far] &&
             m_node_liquid_weight_z[bucket_idx_far][node_idx_far])) {
          m_node_global_pressure_indices[bucket_idx][i] = k++;
        }
      }
    }

    num_effective_nodes[bucket_idx] = k;
  });

  //        std::cout << "[solvePressure: 1]" << std::endl;
  std::partial_sum(num_effective_nodes.begin(), num_effective_nodes.end(),
                   num_effective_nodes.begin());

  const int total_num_nodes =
      num_effective_nodes[num_effective_nodes.size() - 1];
  if (total_num_nodes == 0) return true;

  std::vector<Vec2i> effective_node_indices(total_num_nodes);
  std::vector<Vec3i> dof_ijk(total_num_nodes);
  std::vector<double> result(total_num_nodes);

  result.assign(total_num_nodes, 0.0);
  m_particle_buckets.for_each_bucket(
      [&](int bucket_idx) { m_node_pressure[bucket_idx].setZero(); });

  if ((int)m_pressure_rhs.size() != total_num_nodes) {
    m_pressure_rhs.resize(total_num_nodes);
    m_pressure_matrix.resize(total_num_nodes);
  }

  m_pressure_matrix.zero();
  //        std::cout << "[solvePressure: 2]" << std::endl;
  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    const int num_nodes_p = m_node_liquid_phi[bucket_idx].size();
    const Vec3i handle = m_particle_buckets.bucket_handle(bucket_idx);

    for (int i = 0; i < num_nodes_p; ++i) {
      if (m_node_global_pressure_indices[bucket_idx][i] >= 0) {
        if (bucket_idx > 0) {
          m_node_global_pressure_indices[bucket_idx][i] +=
              num_effective_nodes[bucket_idx - 1];
        }

        const int global_idx = m_node_global_pressure_indices[bucket_idx][i];
        effective_node_indices[global_idx] = Vec2i(bucket_idx, i);

        const Vec3i& local_handle = getNodeHandle(i);
        dof_ijk[global_idx] = Vec3i(handle(0) * m_num_nodes + local_handle(0),
                                    handle(1) * m_num_nodes + local_handle(1),
                                    handle(2) * m_num_nodes + local_handle(2));
      }
    }
  });
  //        std::cout << "[solvePressure: 3]" << std::endl;
  const Scalar dx = getDx();
  const Scalar coeff = m_dt / (dx * dx);
  const Scalar dV = dx * dx * dx;
  const Scalar iD = getInverseDCoeff();

  for_each(0, total_num_nodes, [&](int dof_idx) {
    const Vec2i& dof_loc = effective_node_indices[dof_idx];
    const int bucket_idx = dof_loc[0];
    const int node_idx = dof_loc[1];

    const Scalar center_phi = m_node_liquid_phi[bucket_idx][node_idx];
    m_pressure_rhs[dof_idx] = 0.0;

    const VecXi& bucket_pn = m_node_pressure_neighbors[bucket_idx];

    const int bucket_idx_left = bucket_pn[node_idx * 12 + 0];
    const int node_idx_left = bucket_pn[node_idx * 12 + 1];

    const int bucket_idx_right = bucket_pn[node_idx * 12 + 2];
    const int node_idx_right = bucket_pn[node_idx * 12 + 3];

    const int bucket_idx_bottom = bucket_pn[node_idx * 12 + 4];
    const int node_idx_bottom = bucket_pn[node_idx * 12 + 5];

    const int bucket_idx_top = bucket_pn[node_idx * 12 + 6];
    const int node_idx_top = bucket_pn[node_idx * 12 + 7];

    const int bucket_idx_near = bucket_pn[node_idx * 12 + 8];
    const int node_idx_near = bucket_pn[node_idx * 12 + 9];

    const int bucket_idx_far = bucket_pn[node_idx * 12 + 10];
    const int node_idx_far = bucket_pn[node_idx * 12 + 11];

    Scalar w_center = 0.0;

    if (bucket_idx_left >= 0 && node_idx_left >= 0 &&
        m_bucket_activated[bucket_idx_left]) {
      const Scalar& w = m_node_liquid_weight_x[bucket_idx_left][node_idx_left];

      w_center += w;

      m_pressure_rhs[dof_idx] +=
          (m_node_vel_fluid_plus_x[bucket_idx_left][node_idx_left] * w) / dx;
      if (w > 0.0) {
        const Scalar inv_density =
            (m_node_mass_fluid_x[bucket_idx_left][node_idx_left] < 1e-20)
                ? 1.0
                : (m_node_vol_fluid_x[bucket_idx_left][node_idx_left] /
                   m_node_mass_fluid_x[bucket_idx_left][node_idx_left]);

        Scalar term = w * coeff * inv_density;
        // check_isnan("pleft", term);

        const int bucket_pressure_left =
            m_node_index_pressure_x[bucket_idx_left][node_idx_left * 4 + 0];
        const int node_pressure_left =
            m_node_index_pressure_x[bucket_idx_left][node_idx_left * 4 + 1];

        if (bucket_pressure_left >= 0 && node_pressure_left >= 0 &&
            m_bucket_activated[bucket_pressure_left]) {
          const Scalar left_phi =
              m_node_liquid_phi[bucket_pressure_left][node_pressure_left];

          if (left_phi < 0.0) {
            const int dof_left =
                m_node_global_pressure_indices[bucket_pressure_left]
                                              [node_pressure_left];
            assert(dof_left >= 0);

            if (dof_left >= 0) {
              m_pressure_matrix.add_to_element(dof_idx, dof_idx, term);
              m_pressure_matrix.add_to_element(dof_idx, dof_left, -term);
            }
          } else {
            const Scalar theta = std::max(fraction_inside(center_phi, left_phi),
                                          m_liquid_info.theta_criterion);
            // check_isnan("pleft_theta", theta);
            m_pressure_matrix.add_to_element(dof_idx, dof_idx, term / theta);
          }
        }
      }
    }

    if (bucket_idx_right >= 0 && node_idx_right >= 0 &&
        m_bucket_activated[bucket_idx_right]) {
      const Scalar& w =
          m_node_liquid_weight_x[bucket_idx_right][node_idx_right];

      w_center += w;

      m_pressure_rhs[dof_idx] -=
          (m_node_vel_fluid_plus_x[bucket_idx_right][node_idx_right] * w) / dx;

      if (w > 0.0) {
        const Scalar inv_density =
            (m_node_mass_fluid_x[bucket_idx_right][node_idx_right] < 1e-20)
                ? 1.0
                : (m_node_vol_fluid_x[bucket_idx_right][node_idx_right] /
                   m_node_mass_fluid_x[bucket_idx_right][node_idx_right]);

        Scalar term = w * coeff * inv_density;
        // check_isnan("pright", term);

        const int bucket_pressure_right =
            m_node_index_pressure_x[bucket_idx_right][node_idx_right * 4 + 2];
        const int node_pressure_right =
            m_node_index_pressure_x[bucket_idx_right][node_idx_right * 4 + 3];

        if (bucket_pressure_right >= 0 && node_pressure_right >= 0 &&
            m_bucket_activated[bucket_pressure_right]) {
          const Scalar right_phi =
              m_node_liquid_phi[bucket_pressure_right][node_pressure_right];

          if (right_phi < 0.0) {
            const int dof_right =
                m_node_global_pressure_indices[bucket_pressure_right]
                                              [node_pressure_right];
            assert(dof_right >= 0);

            if (dof_right >= 0) {
              m_pressure_matrix.add_to_element(dof_idx, dof_idx, term);
              m_pressure_matrix.add_to_element(dof_idx, dof_right, -term);
            }
          } else {
            const Scalar theta =
                std::max(fraction_inside(center_phi, right_phi),
                         m_liquid_info.theta_criterion);
            // check_isnan("pright_theta", theta);
            m_pressure_matrix.add_to_element(dof_idx, dof_idx, term / theta);
          }
        }
      }
    }

    if (bucket_idx_bottom >= 0 && node_idx_bottom >= 0 &&
        m_bucket_activated[bucket_idx_bottom]) {
      const Scalar& w =
          m_node_liquid_weight_y[bucket_idx_bottom][node_idx_bottom];

      w_center += w;

      m_pressure_rhs[dof_idx] +=
          (m_node_vel_fluid_plus_y[bucket_idx_bottom][node_idx_bottom] * w) /
          dx;

      if (w > 0.0) {
        const Scalar inv_density =
            (m_node_mass_fluid_y[bucket_idx_bottom][node_idx_bottom] < 1e-20)
                ? 1.0
                : (m_node_vol_fluid_y[bucket_idx_bottom][node_idx_bottom] /
                   m_node_mass_fluid_y[bucket_idx_bottom][node_idx_bottom]);

        Scalar term = w * coeff * inv_density;
        // check_isnan("pbottom", term);

        const int bucket_pressure_bottom =
            m_node_index_pressure_y[bucket_idx_bottom][node_idx_bottom * 4 + 0];
        const int node_pressure_bottom =
            m_node_index_pressure_y[bucket_idx_bottom][node_idx_bottom * 4 + 1];

        if (bucket_pressure_bottom >= 0 && node_pressure_bottom >= 0 &&
            m_bucket_activated[bucket_pressure_bottom]) {
          const Scalar bottom_phi =
              m_node_liquid_phi[bucket_pressure_bottom][node_pressure_bottom];

          if (bottom_phi < 0.0) {
            const int dof_bottom =
                m_node_global_pressure_indices[bucket_pressure_bottom]
                                              [node_pressure_bottom];
            assert(dof_bottom >= 0);

            if (dof_bottom >= 0) {
              m_pressure_matrix.add_to_element(dof_idx, dof_idx, term);
              m_pressure_matrix.add_to_element(dof_idx, dof_bottom, -term);
            }
          } else {
            const Scalar theta =
                std::max(fraction_inside(center_phi, bottom_phi),
                         m_liquid_info.theta_criterion);
            // check_isnan("pbottom_theta", theta);
            m_pressure_matrix.add_to_element(dof_idx, dof_idx, term / theta);
          }
        }
      }
    }

    if (bucket_idx_top >= 0 && node_idx_top >= 0 &&
        m_bucket_activated[bucket_idx_top]) {
      const Scalar& w = m_node_liquid_weight_y[bucket_idx_top][node_idx_top];

      w_center += w;

      m_pressure_rhs[dof_idx] -=
          (m_node_vel_fluid_plus_y[bucket_idx_top][node_idx_top] * w) / dx;

      if (w > 0.0) {
        const Scalar inv_density =
            (m_node_mass_fluid_y[bucket_idx_top][node_idx_top] < 1e-20)
                ? 1.0
                : (m_node_vol_fluid_y[bucket_idx_top][node_idx_top] /
                   m_node_mass_fluid_y[bucket_idx_top][node_idx_top]);

        Scalar term = w * coeff * inv_density;
        // check_isnan("ptop", term);

        const int bucket_pressure_top =
            m_node_index_pressure_y[bucket_idx_top][node_idx_top * 4 + 2];
        const int node_pressure_top =
            m_node_index_pressure_y[bucket_idx_top][node_idx_top * 4 + 3];

        if (bucket_pressure_top >= 0 && node_pressure_top >= 0 &&
            m_bucket_activated[bucket_pressure_top]) {
          const Scalar top_phi =
              m_node_liquid_phi[bucket_pressure_top][node_pressure_top];

          if (top_phi < 0.0) {
            const int dof_top =
                m_node_global_pressure_indices[bucket_pressure_top]
                                              [node_pressure_top];
            assert(dof_top >= 0);

            if (dof_top >= 0) {
              m_pressure_matrix.add_to_element(dof_idx, dof_idx, term);
              m_pressure_matrix.add_to_element(dof_idx, dof_top, -term);
            }
          } else {
            const Scalar theta = std::max(fraction_inside(center_phi, top_phi),
                                          m_liquid_info.theta_criterion);
            // check_isnan("ptop_theta", theta);
            m_pressure_matrix.add_to_element(dof_idx, dof_idx, term / theta);
          }
        }
      }
    }

    if (bucket_idx_near >= 0 && node_idx_near >= 0 &&
        m_bucket_activated[bucket_idx_near]) {
      const Scalar& w = m_node_liquid_weight_z[bucket_idx_near][node_idx_near];

      w_center += w;

      m_pressure_rhs[dof_idx] +=
          (m_node_vel_fluid_plus_z[bucket_idx_near][node_idx_near] * w) / dx;

      if (w > 0.0) {
        const Scalar inv_density =
            (m_node_mass_fluid_z[bucket_idx_near][node_idx_near] < 1e-20)
                ? 1.0
                : (m_node_vol_fluid_z[bucket_idx_near][node_idx_near] /
                   m_node_mass_fluid_z[bucket_idx_near][node_idx_near]);

        Scalar term = w * coeff * inv_density;
        // check_isnan("pnear", term);

        const int bucket_pressure_near =
            m_node_index_pressure_z[bucket_idx_near][node_idx_near * 4 + 0];
        const int node_pressure_near =
            m_node_index_pressure_z[bucket_idx_near][node_idx_near * 4 + 1];

        if (bucket_pressure_near >= 0 && node_pressure_near >= 0 &&
            m_bucket_activated[bucket_pressure_near]) {
          const Scalar near_phi =
              m_node_liquid_phi[bucket_pressure_near][node_pressure_near];

          if (near_phi < 0.0) {
            const int dof_near =
                m_node_global_pressure_indices[bucket_pressure_near]
                                              [node_pressure_near];
            assert(dof_near >= 0);

            if (dof_near >= 0) {
              m_pressure_matrix.add_to_element(dof_idx, dof_idx, term);
              m_pressure_matrix.add_to_element(dof_idx, dof_near, -term);
            }
          } else {
            const Scalar theta = std::max(fraction_inside(center_phi, near_phi),
                                          m_liquid_info.theta_criterion);
            // check_isnan("pnear_theta", theta);
            m_pressure_matrix.add_to_element(dof_idx, dof_idx, term / theta);
          }
        }
      }
    }

    if (bucket_idx_far >= 0 && node_idx_far >= 0 &&
        m_bucket_activated[bucket_idx_far]) {
      const Scalar& w = m_node_liquid_weight_z[bucket_idx_far][node_idx_far];

      w_center += w;

      m_pressure_rhs[dof_idx] -=
          (m_node_vel_fluid_plus_z[bucket_idx_far][node_idx_far] * w) / dx;

      if (w > 0.0) {
        const Scalar inv_density =
            (m_node_mass_fluid_z[bucket_idx_far][node_idx_far] < 1e-20)
                ? 1.0
                : (m_node_vol_fluid_z[bucket_idx_far][node_idx_far] /
                   m_node_mass_fluid_z[bucket_idx_far][node_idx_far]);

        Scalar term = w * coeff * inv_density;
        // check_isnan("pfar", term);

        const int bucket_pressure_far =
            m_node_index_pressure_z[bucket_idx_far][node_idx_far * 4 + 2];
        const int node_pressure_far =
            m_node_index_pressure_z[bucket_idx_far][node_idx_far * 4 + 3];

        if (bucket_pressure_far >= 0 && node_pressure_far >= 0 &&
            m_bucket_activated[bucket_pressure_far]) {
          const Scalar far_phi =
              m_node_liquid_phi[bucket_pressure_far][node_pressure_far];

          if (far_phi < 0.0) {
            const int dof_far =
                m_node_global_pressure_indices[bucket_pressure_far]
                                              [node_pressure_far];
            assert(dof_far >= 0);

            if (dof_far >= 0) {
              m_pressure_matrix.add_to_element(dof_idx, dof_idx, term);
              m_pressure_matrix.add_to_element(dof_idx, dof_far, -term);
            }
          } else {
            const Scalar theta = std::max(fraction_inside(center_phi, far_phi),
                                          m_liquid_info.theta_criterion);
            // check_isnan("pfar_theta", theta);
            m_pressure_matrix.add_to_element(dof_idx, dof_idx, term / theta);
          }
        }
      }
    }

    w_center = std::min(1.0, w_center / 6.0);
    Scalar center_J = m_node_vol_change_p[bucket_idx][node_idx];

    if (center_J > 0.0 && w_center > 0.0) {
      const VecXx& color = m_node_components[bucket_idx].segment(
          node_idx * m_liquid_info.num_components,
          m_liquid_info.num_components);
      const Scalar lambda = getBulkModulus(color) * 0.5;
      const Scalar compliance = 1.0 / (lambda * m_dt);

      center_J = clamp(center_J, 0.99, 1.01);

      const Scalar center_coeff = compliance * w_center / center_J;
      // check_isnan("pcenter", center_coeff);

      m_pressure_matrix.add_to_element(dof_idx, dof_idx, center_coeff);

      const Scalar p_prev = -lambda * (center_J - 1.0 / center_J);

      // check_isnan("p_prev", p_prev);
      m_pressure_rhs[dof_idx] += center_coeff * p_prev;
    }

    if (m_liquid_info.use_varying_volume_fraction && w_center > 0.0) {
      const auto& node_vertex_p = m_node_vertex_p[bucket_idx][node_idx];
      const Vec3x np = getNodePosP(bucket_idx, node_idx);

      for (auto& pair : node_vertex_p) {
        const int pidx = pair.first;

        auto& weights = m_vertex_weights[pidx];
        const int strand_idx = m_strands_global_local[pidx].first;
        const int local_idx = m_strands_global_local[pidx].second;
        const Vec3x pos = m_strands[strand_idx]->getVertex(local_idx);
        const Scalar inv_V_c =
            1.0 / ((1.0 - m_node_elastic_vf_p[bucket_idx][node_idx]) * dV);

        // check_isnan("inv_V_c", inv_V_c);
        // Note: This is the (negative) divergence at grid node! So we use (pos
        // - np) instead of otherwise!
        m_pressure_rhs[dof_idx] -=
            inv_V_c * weights(pair.second, 3) * iD *
            m_elastic_liquid_diffs.segment<3>(pidx * 3).dot(pos - np) *
            w_center;
      }
    }
  });

  double tolerance;
  int iterations;

  bool success =
      AMGPCGSolveSparse(m_pressure_matrix, m_pressure_rhs, result, dof_ijk,
                        1e-8, 1000, tolerance, iterations, ni, nj, nk);
  if (!success) std::cout << "WARNING: Pressure solve failed!" << std::endl;

  std::cout << "[amg pcg total iter: " << iterations << ", res: " << tolerance
            << "]" << std::endl;

  if (iterations >= m_liquid_info.pcg_max_iters) {
    std::cerr << "WARNING: AMG PCG solve failed!" << std::endl;
    std::ofstream ofs("info_amgpcg.m");

    ofs << "rhs=[";
    for (Scalar s : m_pressure_rhs) {
      ofs << s << "; ";
    }
    ofs << "];" << std::endl << std::endl;

    m_pressure_matrix.write_matlab(ofs, "A");
    ofs << std::endl;
    ofs.flush();
    ofs.close();
    exit(-1);
  }

#ifdef CHECK_AMGPCG_RESULT
  std::vector<Scalar> tmp = m_pressure_rhs;

  bridson::multiply_and_subtract(m_pressure_matrix, result, tmp);

  Scalar residual = Eigen::Map<VecXx>((Scalar*)&tmp[0], tmp.size()).norm();

  Scalar len_rhs = Eigen::Map<VecXx>((Scalar*)&rhs[0], rhs.size()).norm();

  std::cout << "[amg pcg check result: " << residual << ", " << len_rhs << ", "
            << (residual / len_rhs) << "]" << std::endl;
#endif
  //        std::cout << "[solvePressure: 8]" << std::endl;
  for_each(0, total_num_nodes, [&](int dof_idx) {
    const Vec2i& dof_loc = effective_node_indices[dof_idx];
    m_node_pressure[dof_loc[0]][dof_loc[1]] = result[dof_idx];
  });
  //        std::cout << "[solvePressure: 9]" << std::endl;
  return true;
}

void LiquidSimulator::extrapolate(vector<VecXuc>& node_valid,
                                  vector<VecXx>& node_vel) {
  const Vec3i offsets[] = {Vec3i(-1, 0, 0), Vec3i(1, 0, 0),  Vec3i(0, -1, 0),
                           Vec3i(0, 1, 0),  Vec3i(0, 0, -1), Vec3i(0, 0, 1)};

  for (int layers = 0; layers < 4; ++layers) {
    std::vector<VecXx> node_vel_copy = node_vel;
    std::vector<VecXuc> node_valid_copy = node_valid;

    m_particle_buckets.for_each_bucket([&](int bucket_idx) {
      if (!m_bucket_activated[bucket_idx]) return;

      const Vec3i bucket_handle = m_particle_buckets.bucket_handle(bucket_idx);

      const VecXuc& bucket_valid = node_valid_copy[bucket_idx];

      for (int k = 0; k < m_num_nodes; ++k)
        for (int j = 0; j < m_num_nodes; ++j)
          for (int i = 0; i < m_num_nodes; ++i) {
            const int center_cell_idx =
                k * m_num_nodes * m_num_nodes + j * m_num_nodes + i;

            if (bucket_valid[center_cell_idx]) continue;

            Scalar sum = 0;
            int count = 0;

            for (int r = 0; r < 6; ++r) {
              Vec3i neigh_idx = Vec3i(i, j, k) + offsets[r];

              Vec3i node_bucket_handle = bucket_handle;
              for (int r = 0; r < 3; ++r) {
                if (neigh_idx(r) < 0) {
                  node_bucket_handle(r)--;
                  neigh_idx(r) += m_num_nodes;
                }
                if (neigh_idx(r) >= m_num_nodes) {
                  node_bucket_handle(r)++;
                  neigh_idx(r) -= m_num_nodes;
                }
              }

              if (node_bucket_handle(0) < 0 ||
                  node_bucket_handle(0) >= m_particle_buckets.dim_size(0) ||
                  node_bucket_handle(1) < 0 ||
                  node_bucket_handle(1) >= m_particle_buckets.dim_size(1) ||
                  node_bucket_handle(2) < 0 ||
                  node_bucket_handle(2) >= m_particle_buckets.dim_size(2))
                continue;

              const int node_bucket_idx =
                  m_particle_buckets.bucket_index(node_bucket_handle);
              if (!m_bucket_activated[node_bucket_idx]) continue;

              const int neigh_cell_idx =
                  neigh_idx(2) * m_num_nodes * m_num_nodes +
                  neigh_idx(1) * m_num_nodes + neigh_idx(0);

              if (!node_valid_copy[node_bucket_idx][neigh_cell_idx]) continue;

              sum += node_vel_copy[node_bucket_idx][neigh_cell_idx];
              ++count;
            }

            if (count > 0) {
              node_vel[bucket_idx][center_cell_idx] = sum / (Scalar)count;
              node_valid[bucket_idx][center_cell_idx] = 1U;
            }
          }
    });
  }
}

bool LiquidSimulator::applyPressure() {
  if (!m_liquid_info.use_implicit_pressure) return true;

  const Scalar dx = getDx();

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    VecXx& bucket_rhs_x = m_node_rhs_fluid_x[bucket_idx];
    VecXx& bucket_rhs_y = m_node_rhs_fluid_y[bucket_idx];
    VecXx& bucket_rhs_z = m_node_rhs_fluid_z[bucket_idx];

    VecXx& bucket_vel_fluid_x = m_node_vel_fluid_plus_x[bucket_idx];
    VecXx& bucket_vel_fluid_y = m_node_vel_fluid_plus_y[bucket_idx];
    VecXx& bucket_vel_fluid_z = m_node_vel_fluid_plus_z[bucket_idx];

    VecXx& bucket_pressure_grad_x = m_node_pressure_grad_x[bucket_idx];
    VecXx& bucket_pressure_grad_y = m_node_pressure_grad_y[bucket_idx];
    VecXx& bucket_pressure_grad_z = m_node_pressure_grad_z[bucket_idx];

    const VecXx& bucket_weight_x = m_node_liquid_weight_x[bucket_idx];
    const VecXx& bucket_weight_y = m_node_liquid_weight_y[bucket_idx];
    const VecXx& bucket_weight_z = m_node_liquid_weight_z[bucket_idx];

    VecXuc& bucket_valid_x = m_node_liquid_valid_x[bucket_idx];
    VecXuc& bucket_valid_y = m_node_liquid_valid_y[bucket_idx];
    VecXuc& bucket_valid_z = m_node_liquid_valid_z[bucket_idx];

    const int num_rhs_x = bucket_rhs_x.size();
    const int num_rhs_y = bucket_rhs_y.size();
    const int num_rhs_z = bucket_rhs_z.size();

    for (int i = 0; i < num_rhs_x; ++i) {
      bucket_valid_x(i) = 0U;
      bucket_pressure_grad_x(i) = 0.0;

      const Scalar vol = m_node_vol_fluid_x[bucket_idx][i];

      if (bucket_weight_x(i) > 0.) {
        const Vec4i& indices =
            m_node_index_pressure_x[bucket_idx].segment<4>(i * 4);

        Scalar left_phi(0), right_phi(0);
        if (indices[0] >= 0 && indices[1] >= 0)
          left_phi = m_node_liquid_phi[indices[0]][indices[1]];

        if (indices[2] >= 0 && indices[3] >= 0)
          right_phi = m_node_liquid_phi[indices[2]][indices[3]];

        if (left_phi < 0.0 || right_phi < 0.0) {
          Scalar theta = 1.;
          if (right_phi >= 0.0 || left_phi >= 0.0)
            theta = std::max(fraction_inside(left_phi, right_phi),
                             m_liquid_info.theta_criterion);

          Scalar pressure_grad = 0;
          if (indices[0] >= 0 && indices[1] >= 0)
            pressure_grad -= m_node_pressure[indices[0]][indices[1]];

          if (indices[2] >= 0 && indices[3] >= 0)
            pressure_grad += m_node_pressure[indices[2]][indices[3]];

          Scalar term = pressure_grad / (dx * theta) * m_dt;

          bucket_rhs_x(i) -= term * vol;
          bucket_pressure_grad_x(i) = pressure_grad / (dx * theta);
          bucket_valid_x(i) = 1U;
        }
      }
    }

    for (int i = 0; i < num_rhs_y; ++i) {
      bucket_valid_y(i) = 0U;
      bucket_pressure_grad_y(i) = 0.0;

      const Scalar vol = m_node_vol_fluid_y[bucket_idx][i];

      if (bucket_weight_y(i) > 0.) {
        const Vec4i& indices =
            m_node_index_pressure_y[bucket_idx].segment<4>(i * 4);

        Scalar bottom_phi(0), top_phi(0);
        if (indices[0] >= 0 && indices[1] >= 0)
          bottom_phi = m_node_liquid_phi[indices[0]][indices[1]];

        if (indices[2] >= 0 && indices[3] >= 0)
          top_phi = m_node_liquid_phi[indices[2]][indices[3]];

        if (bottom_phi < 0.0 || top_phi < 0.0) {
          Scalar theta = 1.;
          if (top_phi >= 0.0 || bottom_phi >= 0.0)
            theta = std::max(fraction_inside(bottom_phi, top_phi),
                             m_liquid_info.theta_criterion);

          Scalar pressure_grad = 0;
          if (indices[0] >= 0 && indices[1] >= 0)
            pressure_grad -= m_node_pressure[indices[0]][indices[1]];

          if (indices[2] >= 0 && indices[3] >= 0)
            pressure_grad += m_node_pressure[indices[2]][indices[3]];

          Scalar term = pressure_grad / (dx * theta) * m_dt;

          bucket_rhs_y(i) -= term * vol;
          bucket_pressure_grad_y(i) = pressure_grad / (dx * theta);
          bucket_valid_y(i) = 1U;
        }
      }
    }

    for (int i = 0; i < num_rhs_z; ++i) {
      bucket_valid_z(i) = 0U;
      bucket_pressure_grad_z(i) = 0.0;

      const Scalar vol = m_node_vol_fluid_z[bucket_idx][i];

      if (bucket_weight_z(i) > 0.) {
        const Vec4i& indices =
            m_node_index_pressure_z[bucket_idx].segment<4>(i * 4);

        Scalar near_phi(0), far_phi(0);
        if (indices[0] >= 0 && indices[1] >= 0)
          near_phi = m_node_liquid_phi[indices[0]][indices[1]];

        if (indices[2] >= 0 && indices[3] >= 0)
          far_phi = m_node_liquid_phi[indices[2]][indices[3]];

        if (near_phi < 0.0 || far_phi < 0.0) {
          Scalar theta = 1.;
          if (far_phi >= 0.0 || near_phi >= 0.0)
            theta = std::max(fraction_inside(near_phi, far_phi),
                             m_liquid_info.theta_criterion);

          Scalar pressure_grad = 0;
          if (indices[0] >= 0 && indices[1] >= 0)
            pressure_grad -= m_node_pressure[indices[0]][indices[1]];

          if (indices[2] >= 0 && indices[3] >= 0)
            pressure_grad += m_node_pressure[indices[2]][indices[3]];

          Scalar term = pressure_grad / (dx * theta) * m_dt;

          bucket_rhs_z(i) -= term * vol;
          bucket_pressure_grad_z(i) = pressure_grad / (dx * theta);
          bucket_valid_z(i) = 1U;
        }
      }
    }
  });

  return true;
}

bool LiquidSimulator::solveViscosity() {
  const int ndof = numParticles() * 3;
  if (!ndof) return true;

  solveForceExplicit(m_node_lhs_fluid_x, m_node_lhs_fluid_y,
                     m_node_lhs_fluid_z);

  if (!m_liquid_info.use_implicit_elasticity || !m_liquid_info.solve_viscosity)
    return true;

  // check_isnan("fluidvel_after_addforce", m_node_vel_fluid_plus_x,
  // m_node_vel_fluid_plus_y, m_node_vel_fluid_plus_z);

  Scalar res_norm_0 = lengthNodeVectors(m_node_rhs_fluid_x, m_node_rhs_fluid_y,
                                        m_node_rhs_fluid_z);

  // check_isnan("res_norm_0 before solveViscosity", res_norm_0);

  if (res_norm_0 > m_liquid_info.shear_pcg_criterion) {
    allocateNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z);
    allocateNodeVectors(m_node_z_x, m_node_z_y, m_node_z_z);
    allocateNodeVectors(m_node_p_x, m_node_p_y, m_node_p_z);
    allocateNodeVectors(m_node_q_x, m_node_q_y, m_node_q_z);

    multiplyHessianMatrix(m_node_vel_fluid_plus_x, m_node_vel_fluid_plus_y,
                          m_node_vel_fluid_plus_z, m_node_r_x, m_node_r_y,
                          m_node_r_z);

    m_particle_buckets.for_each_bucket([&](int bucket_idx) {
      if (!m_bucket_activated[bucket_idx]) return;

      m_node_r_x[bucket_idx] =
          m_node_rhs_fluid_x[bucket_idx] - m_node_r_x[bucket_idx];
      m_node_r_y[bucket_idx] =
          m_node_rhs_fluid_y[bucket_idx] - m_node_r_y[bucket_idx];
      m_node_r_z[bucket_idx] =
          m_node_rhs_fluid_z[bucket_idx] - m_node_r_z[bucket_idx];
    });

    Scalar res_norm =
        lengthNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z) / res_norm_0;

    // check_isnan("res_norm before solveViscosity", res_norm);

    int iter = 0;

    if (res_norm < m_liquid_info.shear_pcg_criterion) {
      std::cout << "[pcg total iter: " << iter << ", res: " << res_norm << "/"
                << m_liquid_info.shear_pcg_criterion
                << ", abs. res: " << (res_norm * res_norm_0) << "/"
                << (m_liquid_info.shear_pcg_criterion * res_norm_0) << "]"
                << std::endl;
    } else {
      solveForceExplicit(m_node_lhs_fluid_x, m_node_lhs_fluid_y,
                         m_node_lhs_fluid_z, m_node_r_x, m_node_r_y, m_node_r_z,
                         m_node_z_x, m_node_z_y, m_node_z_z);

      m_particle_buckets.for_each_bucket([&](int bucket_idx) {
        if (!m_bucket_activated[bucket_idx]) return;

        m_node_p_x[bucket_idx] = m_node_z_x[bucket_idx];
        m_node_p_y[bucket_idx] = m_node_z_y[bucket_idx];
        m_node_p_z[bucket_idx] = m_node_z_z[bucket_idx];
      });

      multiplyHessianMatrix(m_node_p_x, m_node_p_y, m_node_p_z, m_node_q_x,
                            m_node_q_y, m_node_q_z);

      Scalar rho = dotNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z,
                                  m_node_z_x, m_node_z_y, m_node_z_z);

      Scalar alpha = rho / dotNodeVectors(m_node_p_x, m_node_p_y, m_node_p_z,
                                          m_node_q_x, m_node_q_y, m_node_q_z);

      m_particle_buckets.for_each_bucket([&](int bucket_idx) {
        if (!m_bucket_activated[bucket_idx]) return;

        m_node_vel_fluid_plus_x[bucket_idx] += m_node_p_x[bucket_idx] * alpha;
        m_node_r_x[bucket_idx] -= m_node_q_x[bucket_idx] * alpha;
        m_node_vel_fluid_plus_y[bucket_idx] += m_node_p_y[bucket_idx] * alpha;
        m_node_r_y[bucket_idx] -= m_node_q_y[bucket_idx] * alpha;
        m_node_vel_fluid_plus_z[bucket_idx] += m_node_p_z[bucket_idx] * alpha;
        m_node_r_z[bucket_idx] -= m_node_q_z[bucket_idx] * alpha;
      });

      res_norm =
          lengthNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z) / res_norm_0;

      const Scalar rho_criterion =
          (m_liquid_info.shear_pcg_criterion * res_norm_0) *
          (m_liquid_info.shear_pcg_criterion * res_norm_0);

      Scalar rho_old, beta;
      for (;
           iter < m_liquid_info.pcg_max_iters &&
           res_norm > m_liquid_info.shear_pcg_criterion && rho > rho_criterion;
           ++iter) {
        rho_old = rho;

        solveForceExplicit(m_node_lhs_fluid_x, m_node_lhs_fluid_y,
                           m_node_lhs_fluid_z, m_node_r_x, m_node_r_y,
                           m_node_r_z, m_node_z_x, m_node_z_y, m_node_z_z);

        rho = dotNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z, m_node_z_x,
                             m_node_z_y, m_node_z_z);

        beta = rho / rho_old;

        m_particle_buckets.for_each_bucket([&](int bucket_idx) {
          if (!m_bucket_activated[bucket_idx]) return;

          m_node_p_x[bucket_idx] =
              m_node_z_x[bucket_idx] + m_node_p_x[bucket_idx] * beta;
          m_node_p_y[bucket_idx] =
              m_node_z_y[bucket_idx] + m_node_p_y[bucket_idx] * beta;
          m_node_p_z[bucket_idx] =
              m_node_z_z[bucket_idx] + m_node_p_z[bucket_idx] * beta;
        });

        multiplyHessianMatrix(m_node_p_x, m_node_p_y, m_node_p_z, m_node_q_x,
                              m_node_q_y, m_node_q_z);

        alpha = rho / dotNodeVectors(m_node_p_x, m_node_p_y, m_node_p_z,
                                     m_node_q_x, m_node_q_y, m_node_q_z);

        m_particle_buckets.for_each_bucket([&](int bucket_idx) {
          if (!m_bucket_activated[bucket_idx]) return;

          m_node_vel_fluid_plus_x[bucket_idx] += m_node_p_x[bucket_idx] * alpha;
          m_node_r_x[bucket_idx] -= m_node_q_x[bucket_idx] * alpha;
          m_node_vel_fluid_plus_y[bucket_idx] += m_node_p_y[bucket_idx] * alpha;
          m_node_r_y[bucket_idx] -= m_node_q_y[bucket_idx] * alpha;
          m_node_vel_fluid_plus_z[bucket_idx] += m_node_p_z[bucket_idx] * alpha;
          m_node_r_z[bucket_idx] -= m_node_q_z[bucket_idx] * alpha;
        });

        res_norm =
            lengthNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z) / res_norm_0;

        if (m_liquid_info.iteration_print_step > 0 &&
            iter % m_liquid_info.iteration_print_step == 0)
          std::cout << "[pcg total iter: " << iter << ", res: " << res_norm
                    << "/" << m_liquid_info.shear_pcg_criterion
                    << ", abs. res: " << (res_norm * res_norm_0) << "/"
                    << (m_liquid_info.shear_pcg_criterion * res_norm_0)
                    << ", rho: " << (rho / (res_norm_0 * res_norm_0)) << "/"
                    << (rho_criterion / (res_norm_0 * res_norm_0))
                    << ", abs. rho: " << rho << "/" << rho_criterion << "]"
                    << std::endl;
      }

      std::cout << "[pcg total iter: " << iter << ", res: " << res_norm << "/"
                << m_liquid_info.shear_pcg_criterion
                << ", abs. res: " << (res_norm * res_norm_0) << "/"
                << (m_liquid_info.shear_pcg_criterion * res_norm_0)
                << ", rho: " << (rho / (res_norm_0 * res_norm_0)) << "/"
                << (rho_criterion / (res_norm_0 * res_norm_0))
                << ", abs. rho: " << rho << "/" << rho_criterion << "]"
                << std::endl;
    }
  }
  return true;
}
bool LiquidSimulator::solveDrag() { return true; }
bool LiquidSimulator::constrainVelocity() {
  //        std::vector< VecXuc > valid_x = m_node_liquid_valid_x;
  //        std::vector< VecXuc > valid_y = m_node_liquid_valid_y;
  //        std::vector< VecXuc > valid_z = m_node_liquid_valid_z;
  //
  // check_isnan("pressure_grad_before_constrain", m_node_pressure_grad_x,
  // m_node_pressure_grad_y, m_node_pressure_grad_z);

  //        extrapolate(valid_x, m_node_pressure_grad_x);
  //        extrapolate(valid_y, m_node_pressure_grad_y);
  //        extrapolate(valid_z, m_node_pressure_grad_z);
  //
  // check_isnan("pressure_grad_after_extrapolate", m_node_pressure_grad_x,
  // m_node_pressure_grad_y, m_node_pressure_grad_z);
  // check_isnan("vel_fluid_before_constrain", m_node_vel_fluid_plus_x,
  // m_node_vel_fluid_plus_y, m_node_vel_fluid_plus_z);

  extrapolate(m_node_liquid_valid_x, m_node_vel_fluid_plus_x);
  extrapolate(m_node_liquid_valid_y, m_node_vel_fluid_plus_y);
  extrapolate(m_node_liquid_valid_z, m_node_vel_fluid_plus_z);

  // check_isnan("vel_fluid_after_extrapolate", m_node_vel_fluid_plus_x,
  // m_node_vel_fluid_plus_y, m_node_vel_fluid_plus_z);

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    VecXx& node_vel_x = m_node_vel_fluid_plus_x[bucket_idx];
    VecXx& node_vel_y = m_node_vel_fluid_plus_y[bucket_idx];
    VecXx& node_vel_z = m_node_vel_fluid_plus_z[bucket_idx];

    const VecXx& bucket_solid_vel_x = m_node_solid_vel_x[bucket_idx];
    const VecXx& bucket_solid_vel_y = m_node_solid_vel_y[bucket_idx];
    const VecXx& bucket_solid_vel_z = m_node_solid_vel_z[bucket_idx];

    const VecXx& bucket_weight_x = m_node_liquid_weight_x[bucket_idx];
    const VecXx& bucket_weight_y = m_node_liquid_weight_y[bucket_idx];
    const VecXx& bucket_weight_z = m_node_liquid_weight_z[bucket_idx];

    const int num_nodes_x = node_vel_x.size();
    const int num_nodes_y = node_vel_y.size();
    const int num_nodes_z = node_vel_z.size();

    for (int i = 0; i < num_nodes_x; ++i) {
      if (bucket_weight_x[i] == 0.0) {
        node_vel_x[i] = bucket_solid_vel_x[i];
      }
    }

    for (int i = 0; i < num_nodes_y; ++i) {
      if (bucket_weight_y[i] == 0.0) {
        node_vel_y[i] = bucket_solid_vel_y[i];
      }
    }

    for (int i = 0; i < num_nodes_z; ++i) {
      if (bucket_weight_z[i] == 0.0) {
        node_vel_z[i] = bucket_solid_vel_z[i];
      }
    }
  });

  // check_isnan("vel_fluid_after_constrain", m_node_vel_fluid_plus_x,
  // m_node_vel_fluid_plus_y, m_node_vel_fluid_plus_z);

  return true;
}
Vec3x LiquidSimulator::getOrigin() const { return m_grid_mincorner; }
void LiquidSimulator::getDimension(int& ni, int& nj, int& nk) const {
  ni = m_particle_buckets.ni;
  nj = m_particle_buckets.nj;
  nk = m_particle_buckets.nk;
}
Scalar LiquidSimulator::getRadius(int i) const { return m_radius(i); }
Vec3x LiquidSimulator::getParticles(int i) const {
  return m_x.segment<3>(i * 3);
}
bool LiquidSimulator::loadParticles() { return true; }
int LiquidSimulator::numParticles() const { return m_x.size() / 3; }
Vec3x LiquidSimulator::get_future_velocity(const Vec3x& position) const {
  const Scalar dx = getDx();
  Vec3x v;
  v(0) = interpolateValue(position, m_node_vel_fluid_plus_x,
                          m_grid_mincorner + Vec3x(0.0, 0.5, 0.5) * dx, 0.0);
  v(1) = interpolateValue(position, m_node_vel_fluid_plus_y,
                          m_grid_mincorner + Vec3x(0.5, 0.0, 0.5) * dx, 0.0);
  v(2) = interpolateValue(position, m_node_vel_fluid_plus_z,
                          m_grid_mincorner + Vec3x(0.5, 0.5, 0.0) * dx, 0.0);

  return v;
}
Scalar LiquidSimulator::get_mass(const Vec3x& position) const {
  const Scalar dx = getDx();
  Vec3x v;
  v(0) = interpolateValue(position, m_node_mass_fluid_x,
                          m_grid_mincorner + Vec3x(0.0, 0.5, 0.5) * dx, 0.0);
  v(1) = interpolateValue(position, m_node_mass_fluid_y,
                          m_grid_mincorner + Vec3x(0.5, 0.0, 0.5) * dx, 0.0);
  v(2) = interpolateValue(position, m_node_mass_fluid_z,
                          m_grid_mincorner + Vec3x(0.5, 0.5, 0.0) * dx, 0.0);

  return v.sum() / 3.0;
}

Scalar LiquidSimulator::get_volume(const Vec3x& position) const {
  const Scalar dx = getDx();
  Vec3x v;
  v(0) = interpolateValue(position, m_node_vol_fluid_x,
                          m_grid_mincorner + Vec3x(0.0, 0.5, 0.5) * dx, 0.0);
  v(1) = interpolateValue(position, m_node_vol_fluid_y,
                          m_grid_mincorner + Vec3x(0.5, 0.0, 0.5) * dx, 0.0);
  v(2) = interpolateValue(position, m_node_vol_fluid_z,
                          m_grid_mincorner + Vec3x(0.5, 0.5, 0.0) * dx, 0.0);

  return v.sum() / 3.0;
}

Vec3x LiquidSimulator::get_velocity(const Vec3x& position) const {
  const Scalar dx = getDx();
  Vec3x v;
  v(0) = interpolateValue(position, m_node_vel_fluid_x,
                          m_grid_mincorner + Vec3x(0.0, 0.5, 0.5) * dx, 0.0);
  v(1) = interpolateValue(position, m_node_vel_fluid_y,
                          m_grid_mincorner + Vec3x(0.5, 0.0, 0.5) * dx, 0.0);
  v(2) = interpolateValue(position, m_node_vel_fluid_z,
                          m_grid_mincorner + Vec3x(0.5, 0.5, 0.0) * dx, 0.0);

  return v;
}
Vec3x LiquidSimulator::get_solid_velocity(const Vec3x& position) const {
  const Scalar dx = getDx();
  Vec3x v;
  v(0) = interpolateValue(position, m_node_solid_vel_x,
                          m_grid_mincorner + Vec3x(0.0, 0.5, 0.5) * dx, 0.0);
  v(1) = interpolateValue(position, m_node_solid_vel_y,
                          m_grid_mincorner + Vec3x(0.5, 0.0, 0.5) * dx, 0.0);
  v(2) = interpolateValue(position, m_node_solid_vel_z,
                          m_grid_mincorner + Vec3x(0.5, 0.5, 0.0) * dx, 0.0);

  return v;
}

Scalar LiquidSimulator::get_liquid_phi(const Vec3x& position) const {
  const Scalar dx = getDx();
  return interpolateValue(position, m_node_liquid_phi,
                          m_grid_mincorner + Vec3x(0.5, 0.5, 0.5) * dx,
                          3.0 * dx);
}

void LiquidSimulator::storeDragCoeff(int strand_idx, int local_idx,
                                     const Scalar& coeff) {
  const int pidx = m_strands_local_global_base[strand_idx] + local_idx;
  m_strands_drag_coeffs(pidx) = coeff;
}

Vec3x LiquidSimulator::get_velocity(int strand_idx, int local_idx) const {
  const int pidx = m_strands_local_global_base[strand_idx] + local_idx;
  const Mat27x4f& weights = m_vertex_weights[pidx];

  auto& indices_x = m_vertex_nodes_x[pidx];
  auto& indices_y = m_vertex_nodes_y[pidx];
  auto& indices_z = m_vertex_nodes_z[pidx];

  Vec3x vel = Vec3x::Zero();
  for (int i = 0; i < indices_x.rows(); ++i) {
    const int node_bucket_idx = indices_x(i, 0);
    const int node_idx = indices_x(i, 1);
    if (node_bucket_idx == -1 || node_idx == -1 ||
        !m_bucket_activated[node_bucket_idx] || weights(i, 0) == 0.0)
      continue;
    vel(0) += m_node_vel_fluid_x[node_bucket_idx](node_idx) * weights(i, 0);
  }

  for (int i = 0; i < indices_y.rows(); ++i) {
    const int node_bucket_idx = indices_y(i, 0);
    const int node_idx = indices_y(i, 1);
    if (node_bucket_idx == -1 || node_idx == -1 ||
        !m_bucket_activated[node_bucket_idx] || weights(i, 1) == 0.0)
      continue;
    vel(1) += m_node_vel_fluid_y[node_bucket_idx](node_idx) * weights(i, 1);
  }

  for (int i = 0; i < indices_z.rows(); ++i) {
    const int node_bucket_idx = indices_z(i, 0);
    const int node_idx = indices_z(i, 1);
    if (node_bucket_idx == -1 || node_idx == -1 ||
        !m_bucket_activated[node_bucket_idx] || weights(i, 2) == 0.0)
      continue;
    vel(2) += m_node_vel_fluid_z[node_bucket_idx](node_idx) * weights(i, 2);
  }

  return vel;
}

VecXx LiquidSimulator::interpolateValue(const Vec3x& pos,
                                        const std::vector<VecXx>& phi,
                                        const Vec3x& phi_ori,
                                        const VecXx& default_val, int subcells,
                                        int N) const {
  const Scalar dx = getDx();
  Vec3x grid_pos = pos - phi_ori;
  Vec3x base_pos = grid_pos / dx * (Scalar)subcells;
  Vec3i base_idx = Vec3i((int)floor(base_pos(0)), (int)floor(base_pos(1)),
                         (int)floor(base_pos(2)));

  VecXx buf[8];
  const int num_sub_nodes = m_num_nodes * subcells;
  for (int t = 0; t < 2; ++t)
    for (int s = 0; s < 2; ++s)
      for (int r = 0; r < 2; ++r) {
        int local_idx = t * 4 + s * 2 + r;
        Vec3i query_idx = base_idx + Vec3i(r, s, t);
        Vec3i bucket_handle =
            Vec3i(query_idx(0) / num_sub_nodes, query_idx(1) / num_sub_nodes,
                  query_idx(2) / num_sub_nodes);

        if (bucket_handle(0) < 0 || bucket_handle(0) >= m_particle_buckets.ni ||
            bucket_handle(1) < 0 || bucket_handle(1) >= m_particle_buckets.nj ||
            bucket_handle(2) < 0 || bucket_handle(2) >= m_particle_buckets.nk) {
          buf[local_idx] = default_val;
          continue;
        }

        const int bucket_idx = m_particle_buckets.bucket_index(bucket_handle);
        if (!m_bucket_activated[bucket_idx]) {
          buf[local_idx] = default_val;
          continue;
        }

        Vec3i node_handle =
            Vec3i(query_idx(0) - bucket_handle(0) * num_sub_nodes,
                  query_idx(1) - bucket_handle(1) * num_sub_nodes,
                  query_idx(2) - bucket_handle(2) * num_sub_nodes);

        const int node_idx = node_handle(2) * num_sub_nodes * num_sub_nodes +
                             node_handle(1) * num_sub_nodes + node_handle(0);

        buf[local_idx] = phi[bucket_idx].segment(node_idx * N, N);
      }

  Vec3x frac = Vec3x(base_pos(0) - (Scalar)base_idx(0),
                     base_pos(1) - (Scalar)base_idx(1),
                     base_pos(2) - (Scalar)base_idx(2));

  return trilerp(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7],
                 frac[0], frac[1], frac[2]);
}

template <int N>
Eigen::Matrix<Scalar, N, 1> LiquidSimulator::interpolateValue(
    const Vec3x& pos, const std::vector<VecXx>& phi, const Vec3x& phi_ori,
    const Eigen::Matrix<Scalar, N, 1>& default_val, int subcells) const {
  const Scalar dx = getDx();
  Vec3x grid_pos = pos - phi_ori;
  Vec3x base_pos = grid_pos / dx * (Scalar)subcells;
  Vec3i base_idx = Vec3i((int)floor(base_pos(0)), (int)floor(base_pos(1)),
                         (int)floor(base_pos(2)));

  Eigen::Matrix<Scalar, N, 1> buf[8];
  const int num_sub_nodes = m_num_nodes * subcells;
  for (int t = 0; t < 2; ++t)
    for (int s = 0; s < 2; ++s)
      for (int r = 0; r < 2; ++r) {
        int local_idx = t * 4 + s * 2 + r;
        Vec3i query_idx = base_idx + Vec3i(r, s, t);
        Vec3i bucket_handle =
            Vec3i(query_idx(0) / num_sub_nodes, query_idx(1) / num_sub_nodes,
                  query_idx(2) / num_sub_nodes);

        if (bucket_handle(0) < 0 || bucket_handle(0) >= m_particle_buckets.ni ||
            bucket_handle(1) < 0 || bucket_handle(1) >= m_particle_buckets.nj ||
            bucket_handle(2) < 0 || bucket_handle(2) >= m_particle_buckets.nk) {
          buf[local_idx] = default_val;
          continue;
        }

        const int bucket_idx = m_particle_buckets.bucket_index(bucket_handle);
        if (!m_bucket_activated[bucket_idx]) {
          buf[local_idx] = default_val;
          continue;
        }

        Vec3i node_handle =
            Vec3i(query_idx(0) - bucket_handle(0) * num_sub_nodes,
                  query_idx(1) - bucket_handle(1) * num_sub_nodes,
                  query_idx(2) - bucket_handle(2) * num_sub_nodes);

        const int node_idx = node_handle(2) * num_sub_nodes * num_sub_nodes +
                             node_handle(1) * num_sub_nodes + node_handle(0);

        buf[local_idx] = phi[bucket_idx].template segment<N>(node_idx * N);
      }

  Vec3x frac = Vec3x(base_pos(0) - (Scalar)base_idx(0),
                     base_pos(1) - (Scalar)base_idx(1),
                     base_pos(2) - (Scalar)base_idx(2));

  return trilerp(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7],
                 frac[0], frac[1], frac[2]);
}

Scalar LiquidSimulator::interpolateValue(const Vec3x& pos,
                                         const std::vector<VecXx>& phi,
                                         const Vec3x& phi_ori,
                                         const Scalar& default_val,
                                         int subcells) const {
  const Scalar dx = getDx();
  Vec3x grid_pos = pos - phi_ori;
  Vec3x base_pos = grid_pos / dx * (Scalar)subcells;
  Vec3i base_idx = Vec3i((int)floor(base_pos(0)), (int)floor(base_pos(1)),
                         (int)floor(base_pos(2)));

  Scalar buf[8];
  const int num_sub_nodes = m_num_nodes * subcells;
  for (int t = 0; t < 2; ++t)
    for (int s = 0; s < 2; ++s)
      for (int r = 0; r < 2; ++r) {
        int local_idx = t * 4 + s * 2 + r;
        Vec3i query_idx = base_idx + Vec3i(r, s, t);
        Vec3i bucket_handle =
            Vec3i(query_idx(0) / num_sub_nodes, query_idx(1) / num_sub_nodes,
                  query_idx(2) / num_sub_nodes);

        if (bucket_handle(0) < 0 || bucket_handle(0) >= m_particle_buckets.ni ||
            bucket_handle(1) < 0 || bucket_handle(1) >= m_particle_buckets.nj ||
            bucket_handle(2) < 0 || bucket_handle(2) >= m_particle_buckets.nk) {
          buf[local_idx] = default_val;
          continue;
        }

        const int bucket_idx = m_particle_buckets.bucket_index(bucket_handle);
        if (!m_bucket_activated[bucket_idx]) {
          buf[local_idx] = default_val;
          continue;
        }

        Vec3i node_handle =
            Vec3i(query_idx(0) - bucket_handle(0) * num_sub_nodes,
                  query_idx(1) - bucket_handle(1) * num_sub_nodes,
                  query_idx(2) - bucket_handle(2) * num_sub_nodes);

        const int node_idx = node_handle(2) * num_sub_nodes * num_sub_nodes +
                             node_handle(1) * num_sub_nodes + node_handle(0);

        buf[local_idx] = phi[bucket_idx][node_idx];
      }

  Vec3x frac = Vec3x(base_pos(0) - (Scalar)base_idx(0),
                     base_pos(1) - (Scalar)base_idx(1),
                     base_pos(2) - (Scalar)base_idx(2));

  return trilerp(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7],
                 frac[0], frac[1], frac[2]);
}

Vec3x LiquidSimulator::interpolateGradient(const Vec3x& pos,
                                           const std::vector<VecXx>& phi,
                                           const Vec3x& phi_ori,
                                           const Scalar& default_val) const {
  const Scalar dx = getDx();
  Vec3x grid_pos = pos - phi_ori;
  Vec3x base_pos = grid_pos / dx;
  Vec3i base_idx = Vec3i((int)floor(base_pos(0)), (int)floor(base_pos(1)),
                         (int)floor(base_pos(2)));

  Scalar buf[8];
  for (int t = 0; t < 2; ++t)
    for (int s = 0; s < 2; ++s)
      for (int r = 0; r < 2; ++r) {
        int local_idx = t * 4 + s * 2 + r;
        Vec3i query_idx = base_idx + Vec3i(r, s, t);
        Vec3i bucket_handle =
            Vec3i(query_idx(0) / m_num_nodes, query_idx(1) / m_num_nodes,
                  query_idx(2) / m_num_nodes);

        if (bucket_handle(0) < 0 || bucket_handle(0) >= m_particle_buckets.ni ||
            bucket_handle(1) < 0 || bucket_handle(1) >= m_particle_buckets.nj ||
            bucket_handle(2) < 0 || bucket_handle(2) >= m_particle_buckets.nk) {
          buf[local_idx] = default_val;
          continue;
        }

        const int bucket_idx = m_particle_buckets.bucket_index(bucket_handle);
        if (!m_bucket_activated[bucket_idx]) {
          buf[local_idx] = default_val;
          continue;
        }

        Vec3i node_handle =
            Vec3i(query_idx(0) - bucket_handle(0) * m_num_nodes,
                  query_idx(1) - bucket_handle(1) * m_num_nodes,
                  query_idx(2) - bucket_handle(2) * m_num_nodes);

        const int node_idx = node_handle(2) * m_num_nodes * m_num_nodes +
                             node_handle(1) * m_num_nodes + node_handle(0);

        buf[local_idx] = phi[bucket_idx][node_idx];
      }

  Vec3x frac = Vec3x(base_pos(0) - (Scalar)base_idx(0),
                     base_pos(1) - (Scalar)base_idx(1),
                     base_pos(2) - (Scalar)base_idx(2));

  const Scalar ddx00 = (buf[1] - buf[0]);
  const Scalar ddx10 = (buf[3] - buf[2]);
  const Scalar ddx01 = (buf[5] - buf[4]);
  const Scalar ddx11 = (buf[7] - buf[6]);
  const Scalar dv_dx = bilerp(ddx00, ddx10, ddx01, ddx11, frac[1], frac[2]);

  const Scalar ddy00 = (buf[2] - buf[0]);
  const Scalar ddy10 = (buf[3] - buf[1]);
  const Scalar ddy01 = (buf[6] - buf[4]);
  const Scalar ddy11 = (buf[7] - buf[5]);
  const Scalar dv_dy = bilerp(ddy00, ddy10, ddy01, ddy11, frac[0], frac[2]);

  const Scalar ddz00 = (buf[4] - buf[0]);
  const Scalar ddz10 = (buf[5] - buf[1]);
  const Scalar ddz01 = (buf[6] - buf[2]);
  const Scalar ddz11 = (buf[7] - buf[3]);
  const Scalar dv_dz = bilerp(ddz00, ddz10, ddz01, ddz11, frac[0], frac[1]);

  Vec3x n = Vec3x(dv_dx, dv_dy, dv_dz) / dx;

  return n;
}

Scalar LiquidSimulator::interpolateValueAndGradient(
    Vec3x& n, const Vec3x& pos, const std::vector<VecXx>& phi,
    const Vec3x& phi_ori, const Scalar& default_val) const {
  const Scalar dx = getDx();
  Vec3x grid_pos = pos - phi_ori;
  Vec3x base_pos = grid_pos / dx;
  Vec3i base_idx = Vec3i((int)floor(base_pos(0)), (int)floor(base_pos(1)),
                         (int)floor(base_pos(2)));

  Scalar buf[8];
  for (int t = 0; t < 2; ++t)
    for (int s = 0; s < 2; ++s)
      for (int r = 0; r < 2; ++r) {
        int local_idx = t * 4 + s * 2 + r;
        Vec3i query_idx = base_idx + Vec3i(r, s, t);
        Vec3i bucket_handle =
            Vec3i(query_idx(0) / m_num_nodes, query_idx(1) / m_num_nodes,
                  query_idx(2) / m_num_nodes);

        if (bucket_handle(0) < 0 || bucket_handle(0) >= m_particle_buckets.ni ||
            bucket_handle(1) < 0 || bucket_handle(1) >= m_particle_buckets.nj ||
            bucket_handle(2) < 0 || bucket_handle(2) >= m_particle_buckets.nk) {
          buf[local_idx] = default_val;
          continue;
        }

        const int bucket_idx = m_particle_buckets.bucket_index(bucket_handle);
        if (!m_bucket_activated[bucket_idx]) {
          buf[local_idx] = default_val;
          continue;
        }

        Vec3i node_handle =
            Vec3i(query_idx(0) - bucket_handle(0) * m_num_nodes,
                  query_idx(1) - bucket_handle(1) * m_num_nodes,
                  query_idx(2) - bucket_handle(2) * m_num_nodes);

        const int node_idx = node_handle(2) * m_num_nodes * m_num_nodes +
                             node_handle(1) * m_num_nodes + node_handle(0);

        buf[local_idx] = phi[bucket_idx][node_idx];
      }

  Vec3x frac = Vec3x(base_pos(0) - (Scalar)base_idx(0),
                     base_pos(1) - (Scalar)base_idx(1),
                     base_pos(2) - (Scalar)base_idx(2));

  const Scalar ddx00 = (buf[1] - buf[0]);
  const Scalar ddx10 = (buf[3] - buf[2]);
  const Scalar ddx01 = (buf[5] - buf[4]);
  const Scalar ddx11 = (buf[7] - buf[6]);
  const Scalar dv_dx = bilerp(ddx00, ddx10, ddx01, ddx11, frac[1], frac[2]);

  const Scalar ddy00 = (buf[2] - buf[0]);
  const Scalar ddy10 = (buf[3] - buf[1]);
  const Scalar ddy01 = (buf[6] - buf[4]);
  const Scalar ddy11 = (buf[7] - buf[5]);
  const Scalar dv_dy = bilerp(ddy00, ddy10, ddy01, ddy11, frac[0], frac[2]);

  const Scalar ddz00 = (buf[4] - buf[0]);
  const Scalar ddz10 = (buf[5] - buf[1]);
  const Scalar ddz01 = (buf[6] - buf[2]);
  const Scalar ddz11 = (buf[7] - buf[3]);
  const Scalar dv_dz = bilerp(ddz00, ddz10, ddz01, ddz11, frac[0], frac[1]);

  n = Vec3x(dv_dx, dv_dy, dv_dz) / dx;

  return trilerp(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7],
                 frac[0], frac[1], frac[2]);
}

void LiquidSimulator::makeParticlesComponentConsistent() {
  const int num_fluid = numParticles();
  for_each(0, num_fluid, [&](int pidx) {
    const Scalar sum = m_components
                           .segment(pidx * m_liquid_info.num_components,
                                    m_liquid_info.num_components)
                           .sum();
    if (sum == 0.0) {
      m_components
          .segment(pidx * m_liquid_info.num_components,
                   m_liquid_info.num_components)
          .setZero();
      m_components(pidx * m_liquid_info.num_components) = 1.0;
    } else {
      m_components.segment(pidx * m_liquid_info.num_components,
                           m_liquid_info.num_components) /= sum;
    }
  });
}

void LiquidSimulator::makeParticlesConsistent() {
  const int num_fluid = numParticles();
  for_each(0, num_fluid, [&](int pidx) {
    const Scalar vol = m_vol(pidx);
    const Scalar r = pow(vol * 0.75 / M_PI, 1.0 / 3.0);
    const VecXx& color = m_components.segment(
        pidx * m_liquid_info.num_components, m_liquid_info.num_components);
    const Scalar m = vol * getDensity(color);

    m_radius(pidx) = r;
    m_m.segment<3>(pidx * 3).setConstant(m);
  });
}

void LiquidSimulator::absorbLiquidFromMesh() {
  const Scalar dx = getDx();

  const int num_fields = m_fields.size();
  const int num_fluid = numParticles();
  // prepare sorter
  m_particle_cells_single.sort(
      num_fluid, [&](int pidx, int& i, int& j, int& k) {
        Vec3x local_x = (m_x.segment<3>(pidx * 3) - m_grid_mincorner) / dx;
        i = (int)floor(local_x(0));
        j = (int)floor(local_x(1));
        k = (int)floor(local_x(2));
      });

  m_mesh_vertex_buckets.for_each_bucket_colored(
      [&](int bucket_idx) {
        const int total_num_cells = m_num_nodes * m_num_nodes * m_num_nodes;

        VecXx reservoir(total_num_cells);
        reservoir.setZero();

        std::vector<int> buf_pidx;

        m_mesh_vertex_buckets.get_bucket(bucket_idx, [&](int pidx) {
          const int mesh_idx = m_mesh_global_local[pidx].first;
          const int local_idx = m_mesh_global_local[pidx].second;

          if (!m_fields[mesh_idx].has_surf_flow) return;

          std::shared_ptr<TriangularMeshFlow> flow =
              m_fields[mesh_idx].mesh_controller->getMeshFlow();

          if (!flow) return;

          std::shared_ptr<TriangularMesh> mesh =
              m_fields[mesh_idx].mesh_controller->getCurrentMesh();

          Vec3x pos = mesh->getVertex(local_idx);

          Scalar phi = get_dense_liquid_signed_distance(pos);

          if (phi >= -1.5 * dx) return;

          Scalar vol = flow->getVertexVol(local_idx);

          // find pos in reservoir
          Vec3x base_x = getNodePos(bucket_idx, 0, 4);
          Vec3x local_x = (mesh->getVertex(local_idx) - base_x) / dx;

          Vec3i ilocal_x = Vec3i(clamp((int)local_x(0), 0, m_num_nodes),
                                 clamp((int)local_x(1), 0, m_num_nodes),
                                 clamp((int)local_x(2), 0, m_num_nodes));
          const int node_idx = ilocal_x(2) * m_num_nodes * m_num_nodes +
                               ilocal_x(1) * m_num_nodes + ilocal_x(0);

          reservoir(node_idx) += vol;
          buf_pidx.push_back(pidx);
        });

        bool has_at_least_one_cell_available = false;

        for (int i = 0; i < total_num_cells; ++i) {
          if (reservoir(i) > 0.0) {
            has_at_least_one_cell_available = true;
          }
        }

        if (!has_at_least_one_cell_available) return;

        Vec3i bucket_handle = m_mesh_vertex_buckets.bucket_handle(bucket_idx);
        Vec3i cell_handle_base = Vec3i(bucket_handle(0) * m_num_nodes,
                                       bucket_handle(1) * m_num_nodes,
                                       bucket_handle(2) * m_num_nodes);

        VecXx old_reservoir = reservoir;
        // for each reservoir cell, firstly try put liquid onto particles
        for (int k = 0; k < m_num_nodes; ++k)
          for (int j = 0; j < m_num_nodes; ++j)
            for (int i = 0; i < m_num_nodes; ++i) {
              int local_idx =
                  k * m_num_nodes * m_num_nodes + j * m_num_nodes + i;
              if (reservoir(local_idx) == 0.0) continue;

              Vec3i cell_handle = cell_handle_base + Vec3i(i, j, k);
              if (!m_particle_cells_single.has_bucket(cell_handle)) continue;

              int cell_idx = m_particle_cells_single.bucket_index(cell_handle);
              int num_parts =
                  m_particle_cells_single.get_neighbor_bucket_particles_size(
                      cell_idx);
              if (!num_parts) continue;

              Scalar vol_each_p = reservoir(local_idx) / (Scalar)num_parts;

              m_particle_cells_single.loop_neighbor_bucket_particles(
                  cell_idx, [&](int pidx, int) {
                    m_vol(pidx) *=
                        (m_rest_vol(pidx) + vol_each_p) / m_rest_vol(pidx);
                    m_rest_vol(pidx) += vol_each_p;

                    return false;
                  });

              reservoir(local_idx) = 0.0;
            }

        // adjust height field
        for (int pidx : buf_pidx) {
          const int mesh_idx = m_mesh_global_local[pidx].first;
          const int local_idx = m_mesh_global_local[pidx].second;

          if (!m_fields[mesh_idx].has_surf_flow) return;

          std::shared_ptr<TriangularMeshFlow> flow =
              m_fields[mesh_idx].mesh_controller->getMeshFlow();

          if (!flow) return;

          Scalar height = flow->getVertexHeight(local_idx);

          std::shared_ptr<TriangularMesh> mesh =
              m_fields[mesh_idx].mesh_controller->getCurrentMesh();

          // find pos in reservoir
          Vec3x base_x = getNodePos(bucket_idx, 0, 4);
          Vec3x local_x = (mesh->getVertex(local_idx) - base_x) / dx;

          Vec3i ilocal_x = Vec3i(clamp((int)local_x(0), 0, m_num_nodes),
                                 clamp((int)local_x(1), 0, m_num_nodes),
                                 clamp((int)local_x(2), 0, m_num_nodes));
          const int node_idx = ilocal_x(2) * m_num_nodes * m_num_nodes +
                               ilocal_x(1) * m_num_nodes + ilocal_x(0);

          if (old_reservoir(node_idx) == 0.0) return;

          height *= reservoir(node_idx) / old_reservoir(node_idx);
          flow->setVertexHeight(local_idx, height);
        }
      },
      3);

  makeParticlesConsistent();
}

void LiquidSimulator::transferBulkLiquidMesh() {
  const int num_fields = m_fields.size();
  if (!num_fields) return;

  bool has_flow = false;
  for (int i = 0; i < num_fields; ++i) {
    has_flow = has_flow || m_fields[i].has_surf_flow;
    if (has_flow) break;
  }

  if (!has_flow) return;

  m_mesh_vertex_buckets.sort(
      m_mesh_global_local.size(), [&](int pidx, int& i, int& j, int& k) {
        const pair<int, int>& local_idx = m_mesh_global_local[pidx];
        const std::shared_ptr<TriangularMesh> mesh =
            m_fields[local_idx.first].mesh_controller->getCurrentMesh();
        const Vec3x pos = mesh->getVertex(local_idx.second);

        i = (int)std::floor((pos(0) - m_grid_mincorner(0)) / m_bucket_size);
        j = (int)std::floor((pos(1) - m_grid_mincorner(1)) / m_bucket_size);
        k = (int)std::floor((pos(2) - m_grid_mincorner(2)) / m_bucket_size);
      });

  const int num_fluid = numParticles();

  if (m_liquid_info.use_liquid_capture && num_fluid) captureBulkLiquidMesh();

  dripFromMesh();
  //
  //        if(num_fluid)
  //            absorbLiquidFromMesh();
}

void LiquidSimulator::captureBulkLiquidMesh() {
  const Scalar h_t = m_liquid_info.typical_flow_thickness * 2.0;
  const Scalar dx = getDx();

  const int num_fluid = numParticles();
  const int num_fields = m_fields.size();
  // prepare flow
  for_each(0, num_fields, [&](int field_idx) {
    std::shared_ptr<TriangularMeshFlow> flow =
        m_fields[field_idx].mesh_controller->getMeshFlow();
    if (!flow) return;

    flow->check_passing_vertices();
  });

  // prepare sorter
  m_particle_cells_single.sort(
      num_fluid, [&](int pidx, int& i, int& j, int& k) {
        Vec3x local_x = (m_x.segment<3>(pidx * 3) - m_grid_mincorner) / dx;
        i = (int)floor(local_x(0));
        j = (int)floor(local_x(1));
        k = (int)floor(local_x(2));
      });

  m_mesh_vertex_buckets.for_each_bucket_particles_colored([&](int pidx, int) {
    const int mesh_idx = m_mesh_global_local[pidx].first;
    const int local_idx = m_mesh_global_local[pidx].second;

    if (!m_fields[mesh_idx].has_surf_flow) return;

    std::shared_ptr<TriangularMeshFlow> flow =
        m_fields[mesh_idx].mesh_controller->getMeshFlow();

    if (!flow || !flow->isVertexPassing(local_idx)) return;

    Scalar height = flow->getVertexHeight(local_idx);

    if (height >= h_t) return;

    std::shared_ptr<TriangularMesh> mesh =
        m_fields[mesh_idx].mesh_controller->getCurrentMesh();

    Scalar area = mesh->getVertexArea(local_idx);
    Scalar vol_to_capture = (h_t - height) * area;

    const Vec3x pos = mesh->getVertex(local_idx);
    const Vec3i ipos = Vec3i((pos(0) - m_grid_mincorner(0)) / dx,
                             (pos(1) - m_grid_mincorner(1)) / dx,
                             (pos(2) - m_grid_mincorner(2)) / dx);

    if (ipos(0) < 0 || ipos(0) >= m_particle_cells_single.ni || ipos(1) < 0 ||
        ipos(1) >= m_particle_cells_single.nj || ipos(2) < 0 ||
        ipos(2) >= m_particle_cells_single.nk)
      return;  // shouldn't reach here!

    const int cell_idx = m_particle_cells_single.bucket_index(ipos);
    const Vec3x u_s = mesh->getDisplacement(local_idx) / m_dt;

    Scalar vol_remain_to_cap = vol_to_capture;

    Vec3x momentum = Vec3x::Zero();
    VecXx color = flow->getCurrentComponents().segment(
                      local_idx * m_liquid_info.num_components,
                      m_liquid_info.num_components) *
                  height * area;

    m_particle_cells_single.loop_neighbor_bucket_particles(
        cell_idx,
        [&](int npidx, int) {
          const Vec3x& np = m_x.segment<3>(npidx * 3);
          Vec3x diff = (pos - np) / dx;
          const Scalar w = N_kernel<2>(diff);
          Scalar vol_capture =
              std::min(m_rest_vol(npidx) * w, vol_remain_to_cap);

          vol_remain_to_cap -= vol_capture;
          m_vol(npidx) *= (m_rest_vol(npidx) - vol_capture) / m_rest_vol(npidx);
          m_rest_vol(npidx) -= vol_capture;

          momentum += m_v.segment<3>(npidx * 3) * vol_capture;
          color += m_components.segment(npidx * m_liquid_info.num_components,
                                        m_liquid_info.num_components) *
                   vol_capture;

          return vol_remain_to_cap <= 0.0;
        },
        2);

    Scalar actual_captured = vol_to_capture - vol_remain_to_cap;

    if (actual_captured <= 0.0) return;

    Scalar actual_height = actual_captured / area + height;

    Vec3x mom_rel = momentum - u_s * actual_captured;

    const std::vector<int>& faces = mesh->getNeighborFacesAtVertex(local_idx);

    if (!faces.size()) return;

    const Vec3x mom_disp = mom_rel / (Scalar)faces.size();
    const Scalar vol_disp = actual_captured / (Scalar)faces.size();

    for (int fidx : faces) {
      Scalar ori_vol = mesh->m_face_areas(fidx) * flow->getFaceHeight(fidx);
      Vec3x ori_mom = flow->getFaceVelocity(fidx) * ori_vol;
      Vec3x new_vel_edge = (ori_mom + mom_disp) / (ori_vol + vol_disp);

      Vec3x normal = mesh->getFaceNormal(fidx);
      new_vel_edge =
          (Mat3x::Identity() - normal * normal.transpose()) * new_vel_edge;
      flow->setFaceVelocity(fidx, new_vel_edge);
    }

    flow->setVertexHeight(local_idx, actual_height);

    make_gibbs_simplex(color);

    flow->setVertexComponents(local_idx, color);
  });

  // remove empty particles
  removeEmptyParticles();
  makeParticlesConsistent();
}

void LiquidSimulator::dripFromMesh() {
  // for each vertex check if its height has reached critial threshold, if so
  // add the extra volume into volume that can be released principle: try absorb
  // extra liquid as much as possible with current particles
  const Scalar dx = getDx();

  const int num_fields = m_fields.size();
  const int num_fluid = numParticles();
  // prepare sorter
  m_particle_cells_single.sort(
      num_fluid, [&](int pidx, int& i, int& j, int& k) {
        Vec3x local_x = (m_x.segment<3>(pidx * 3) - m_grid_mincorner) / dx;
        i = (int)floor(local_x(0));
        j = (int)floor(local_x(1));
        k = (int)floor(local_x(2));
      });

  const Scalar rad =
      defaultRadiusMultiplier() * dx * m_liquid_info.particle_cell_multiplier;
  const Scalar pvol = 4.0 / 3.0 * M_PI * rad * rad * rad;

  struct Sample {
    Vec3x p;
    Vec3x v;
    VecXx c;
  };

  vector<vector<Sample> > bucket_samples(m_mesh_vertex_buckets.ni *
                                         m_mesh_vertex_buckets.nj *
                                         m_mesh_vertex_buckets.nk);

  //        std::cout << getTotalVolMeshFlows() << " + " <<
  //        getTotalVolParticles() << std::endl;

  m_mesh_vertex_buckets.for_each_bucket_colored(
      [&](int bucket_idx) {
        const int total_num_cells = m_num_nodes * m_num_nodes * m_num_nodes;

        VecXx reservoir(total_num_cells);
        reservoir.setZero();

        VecXx average_pos(total_num_cells * 3);
        average_pos.setZero();

        VecXx average_vel(total_num_cells * 3);
        average_vel.setZero();

        VecXx average_color(total_num_cells * m_liquid_info.num_components);
        average_color.setZero();

        m_mesh_vertex_buckets.get_bucket(bucket_idx, [&](int pidx) {
          const int mesh_idx = m_mesh_global_local[pidx].first;
          const int local_idx = m_mesh_global_local[pidx].second;

          if (!m_fields[mesh_idx].has_surf_flow) return;

          std::shared_ptr<TriangularMeshFlow> flow =
              m_fields[mesh_idx].mesh_controller->getMeshFlow();

          if (!flow) return;

          Scalar height = flow->getVertexHeight(local_idx);

          if (height <= m_liquid_info.mesh_flow_critical_height) return;

          std::shared_ptr<TriangularMesh> mesh =
              m_fields[mesh_idx].mesh_controller->getCurrentMesh();
          Scalar area = mesh->getVertexArea(local_idx);

          Scalar extra_vol =
              (height - m_liquid_info.mesh_flow_critical_height) * area;

          // find pos in reservoir
          Vec3x base_x = getNodePos(bucket_idx, 0, 4);
          Vec3x local_x = (mesh->getVertex(local_idx) - base_x) / dx;

          Vec3i ilocal_x = Vec3i(clamp((int)local_x(0), 0, m_num_nodes),
                                 clamp((int)local_x(1), 0, m_num_nodes),
                                 clamp((int)local_x(2), 0, m_num_nodes));
          const int node_idx = ilocal_x(2) * m_num_nodes * m_num_nodes +
                               ilocal_x(1) * m_num_nodes + ilocal_x(0);

          reservoir(node_idx) += extra_vol;
          average_pos.segment<3>(node_idx * 3) +=
              mesh->getVertex(local_idx) * extra_vol;
          average_vel.segment<3>(node_idx * 3) +=
              flow->getVertexVelocity(local_idx) * extra_vol;
          average_color.segment(node_idx * m_liquid_info.num_components,
                                m_liquid_info.num_components) +=
              flow->getCurrentComponents().segment(
                  local_idx * m_liquid_info.num_components,
                  m_liquid_info.num_components) *
              extra_vol;
        });

        bool has_at_least_one_cell_available = false;

        for (int i = 0; i < total_num_cells; ++i) {
          if (reservoir(i) > 0.0) {
            has_at_least_one_cell_available = true;
            average_pos.segment<3>(i * 3) /= reservoir(i);
            average_vel.segment<3>(i * 3) /= reservoir(i);

            make_gibbs_simplex(average_color);
          }
        }

        if (!has_at_least_one_cell_available) return;

        Vec3i bucket_handle = m_mesh_vertex_buckets.bucket_handle(bucket_idx);
        Vec3i cell_handle_base = Vec3i(bucket_handle(0) * m_num_nodes,
                                       bucket_handle(1) * m_num_nodes,
                                       bucket_handle(2) * m_num_nodes);

        VecXx old_reservoir = reservoir;

        // for remaining reservoir cells, check if we can generate particles
        for (int k = 0; k < m_num_nodes; ++k)
          for (int j = 0; j < m_num_nodes; ++j)
            for (int i = 0; i < m_num_nodes; ++i) {
              int local_idx =
                  k * m_num_nodes * m_num_nodes + j * m_num_nodes + i;
              if (reservoir(local_idx) == 0.0 || reservoir(local_idx) < pvol)
                continue;

              const int num_releases = (int)(reservoir(local_idx) / pvol);
              const Scalar actual_released_vol = num_releases * pvol;

              reservoir(local_idx) -= actual_released_vol;

              const Scalar radius_released_vol =
                  pow(actual_released_vol * 0.75 / M_PI, 1.0 / 3.0) * 2.0;

              // generate particles
              Mat3x M = Mat3x::Random();
              Mat3x Q, R;
              QRDecompose<Scalar, 3>(M, Q, R);

              const Vec3x pos = average_pos.segment<3>(local_idx * 3);

              for (int i = 0; i < num_releases; ++i) {
                Sample s;
                Vec3x actual_pos =
                    pos + Q * m_sphere_pattern[num_releases].segment<3>(i * 3) *
                              radius_released_vol;

                // pos is available
                s.v = average_vel.segment<3>(local_idx * 3);
                s.p = actual_pos;
                s.c = average_color;
                bucket_samples[bucket_idx].push_back(s);
              }
            }

        // adjust height field
        m_mesh_vertex_buckets.get_bucket(bucket_idx, [&](int pidx) {
          const int mesh_idx = m_mesh_global_local[pidx].first;
          const int local_idx = m_mesh_global_local[pidx].second;

          if (!m_fields[mesh_idx].has_surf_flow) return;

          std::shared_ptr<TriangularMeshFlow> flow =
              m_fields[mesh_idx].mesh_controller->getMeshFlow();

          if (!flow) return;

          Scalar height = flow->getVertexHeight(local_idx);

          if (height <= m_liquid_info.mesh_flow_critical_height) return;

          std::shared_ptr<TriangularMesh> mesh =
              m_fields[mesh_idx].mesh_controller->getCurrentMesh();

          // find pos in reservoir
          Vec3x base_x = getNodePos(bucket_idx, 0, 4);
          Vec3x local_x = (mesh->getVertex(local_idx) - base_x) / dx;

          Vec3i ilocal_x = Vec3i(clamp((int)local_x(0), 0, m_num_nodes),
                                 clamp((int)local_x(1), 0, m_num_nodes),
                                 clamp((int)local_x(2), 0, m_num_nodes));
          const int node_idx = ilocal_x(2) * m_num_nodes * m_num_nodes +
                               ilocal_x(1) * m_num_nodes + ilocal_x(0);

          if (old_reservoir(node_idx) == 0.0) return;

          height = m_liquid_info.mesh_flow_critical_height +
                   (height - m_liquid_info.mesh_flow_critical_height) *
                       reservoir(node_idx) / old_reservoir(node_idx);
          flow->setVertexHeight(local_idx, height);
        });
      },
      3);

  vector<Sample> new_samples;
  parallel_concatenate(new_samples, bucket_samples);

  const int num_new_samples = new_samples.size();
  if (!num_new_samples) return;

  {
    LockGuard lock(m_particle_mutex);

    const int num_old_parts = numParticles();
    conservativeResizeParticles(num_old_parts + num_new_samples);

    for_each(0, num_new_samples, [&](int i) {
      const int part_idx = num_old_parts + i;
      m_x.segment<3>(part_idx * 3) = new_samples[i].p;
      m_v.segment<3>(part_idx * 3) = new_samples[i].v;
      m_m.segment<3>(part_idx * 3)
          .setConstant(pvol * getDensity(new_samples[i].c));
      m_vol(part_idx) = m_rest_vol(part_idx) = pvol;
      m_radius(part_idx) = rad;
      m_particle_group[part_idx] = 0;
      m_B.block<3, 3>(part_idx * 3, 0).setZero();
      m_Fe.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_b.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_b_trial.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_b_plus.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_Fe_plus.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_proj_func.segment<3>(part_idx * 3) = Vec3x(1, 0, 0);
      m_J(part_idx) = 1.0;
      m_weakened[part_idx] = false;
      m_components.segment(part_idx * m_liquid_info.num_components,
                           m_liquid_info.num_components) = new_samples[i].c;
    });
  }
  makeParticlesConsistent();

  //        std::cout << getTotalVolMeshFlows() << " + " <<
  //        getTotalVolParticles() << std::endl;
}

Vec3i LiquidSimulator::getGridDimensions() const {
  return Vec3i(m_particle_buckets.ni * m_num_nodes,
               m_particle_buckets.nj * m_num_nodes,
               m_particle_buckets.nk * m_num_nodes);
}

void LiquidSimulator::captureBulkLiquid() {
  if (!m_liquid_info.use_liquid_capture) return;

  const Scalar h_t = m_liquid_info.typical_flow_thickness;
  const Scalar dx = getDx();

  const int num_fluid = numParticles();

  if (num_fluid == 0) return;

  const int num_strand = m_strands.size();
  if (!num_strand) return;

  // prepare strands
  for_each(0, num_strand, [&](int strand_idx) {
    ElasticStrand* strand = m_strands[strand_idx];
    strand->dynamics().checkPassingVertices();
  });

  // prepare sorter (possibly terminated by previous process)
  m_particle_cells_single.sort(
      num_fluid, [&](int pidx, int& i, int& j, int& k) {
        Vec3x local_x = (m_x.segment<3>(pidx * 3) - m_grid_mincorner) / dx;
        i = (int)floor(local_x(0));
        j = (int)floor(local_x(1));
        k = (int)floor(local_x(2));
      });

  // update vertex bucket (since vertex pos has been updated!)
  m_vertex_buckets.sort(
      m_strands_global_local.size(), [&](int pidx, int& i, int& j, int& k) {
        const pair<int, int>& local_idx = m_strands_global_local[pidx];
        const ElasticStrand* strand = m_strands[local_idx.first];
        const Vec3x pos = strand->getVertex(local_idx.second);
        i = (int)std::floor((pos(0) - m_grid_mincorner(0)) / m_bucket_size);
        j = (int)std::floor((pos(1) - m_grid_mincorner(1)) / m_bucket_size);
        k = (int)std::floor((pos(2) - m_grid_mincorner(2)) / m_bucket_size);
      });

  // for each passing vertex, find neighbor particles and capture
  m_vertex_buckets.for_each_bucket_particles_colored([&](int pidx, int) {
    const int strand_idx = m_strands_global_local[pidx].first;
    const int local_idx = m_strands_global_local[pidx].second;

    // check if near to liquid
    const StrandDynamicTraits& dyn = m_strands[strand_idx]->dynamics();

    if (!dyn.isVertexPassing(local_idx)) return;

    const Scalar len = m_strands[strand_idx]->getVoronoiLength(local_idx);
    const Scalar area = m_strands[strand_idx]->getCurrentFlowDOFArea(local_idx);

    const Scalar ra = m_strands[strand_idx]->getRadiusA(local_idx);
    const Scalar rb = m_strands[strand_idx]->getRadiusB(local_idx);

    // determain volume to be transferred, ignore the ones that has more flow
    // than we should capture
    const Scalar capture_area = M_PI * h_t * (h_t + ra + rb);
    if (area >= capture_area) return;

    const Scalar vol = (capture_area - area) * len;

    const Vec3x pos = m_strands[strand_idx]->getVertex(local_idx);
    const Vec3i ipos = Vec3i((pos(0) - m_grid_mincorner(0)) / dx,
                             (pos(1) - m_grid_mincorner(1)) / dx,
                             (pos(2) - m_grid_mincorner(2)) / dx);

    if (ipos(0) < 0 || ipos(0) >= m_particle_cells_single.ni || ipos(1) < 0 ||
        ipos(1) >= m_particle_cells_single.nj || ipos(2) < 0 ||
        ipos(2) >= m_particle_cells_single.nk)
      return;  // shouldn't reach here!

    const int cell_idx = m_particle_cells_single.bucket_index(ipos);

    Scalar vol_remain_to_cap = vol;
    const Vec3x u_s = dyn.getDisplacement(local_idx) / m_dt;

    Vec3x momentum = Vec3x::Zero();
    VecXx color = m_strands[strand_idx]->getStepper()->flowComponents().segment(
                      local_idx * m_liquid_info.num_components,
                      m_liquid_info.num_components) *
                  area * len;

    m_particle_cells_single.loop_neighbor_bucket_particles(
        cell_idx,
        [&](int npidx, int ncidx) {
          const Vec3x& np = m_x.segment<3>(npidx * 3);
          Vec3x diff = (pos - np) / dx;
          const Scalar w = N_kernel<2>(diff);
          Scalar vol_capture =
              std::min(m_rest_vol(npidx) * w, vol_remain_to_cap);

          vol_remain_to_cap -= vol_capture;
          m_vol(npidx) *= (m_rest_vol(npidx) - vol_capture) / m_rest_vol(npidx);
          m_rest_vol(npidx) -= vol_capture;

          momentum += m_v.segment<3>(npidx * 3) * vol_capture;
          color += m_components.segment(npidx * m_liquid_info.num_components,
                                        m_liquid_info.num_components) *
                   vol_capture;

          return vol_remain_to_cap <= 0.0;
        },
        2);

    Scalar actual_captured = vol - vol_remain_to_cap;

    if (actual_captured == 0.0) return;

    Scalar actual_area = actual_captured / len + area;

    m_strands[strand_idx]->setCurrentFlowDOFArea(local_idx, actual_area);

    make_gibbs_simplex(color);

    m_strands[strand_idx]->getStepper()->flowComponents().segment(
        local_idx * m_liquid_info.num_components,
        m_liquid_info.num_components) = color;

    Vec3x mom_rel = momentum - u_s * actual_captured;

    const int num_verts = m_strands[strand_idx]->getNumVertices();

    if (local_idx != 0) {
      Scalar ori_vol = m_strands[strand_idx]->getCurrentSurfaceFlowVolumeAtEdge(
          local_idx - 1);
      Scalar mom_edge = mom_rel.dot(
          m_strands[strand_idx]->getEdgeVector(local_idx - 1).normalized());
      Scalar ori_mom_edge =
          m_strands[strand_idx]->getStepper()->getCurrentFlowVelocity()(
              local_idx - 1) *
          ori_vol;

      Scalar fraction = (local_idx == num_verts - 1) ? 1.0 : 0.5;
      Scalar new_vel_edge = (ori_mom_edge + mom_edge * fraction) /
                            (ori_vol + actual_captured * fraction);

      m_strands[strand_idx]->getStepper()->setCurrentFlowVelocity(local_idx - 1,
                                                                  new_vel_edge);
    }

    if (local_idx != num_verts - 1) {
      Scalar ori_vol =
          m_strands[strand_idx]->getCurrentSurfaceFlowVolumeAtEdge(local_idx);
      Scalar mom_edge = mom_rel.dot(
          m_strands[strand_idx]->getEdgeVector(local_idx).normalized());
      Scalar ori_mom_edge =
          m_strands[strand_idx]->getStepper()->getCurrentFlowVelocity()(
              local_idx) *
          ori_vol;

      Scalar fraction = (local_idx == 0) ? 1.0 : 0.5;
      Scalar new_vel_edge = (ori_mom_edge + mom_edge * fraction) /
                            (ori_vol + actual_captured * fraction);

      m_strands[strand_idx]->getStepper()->setCurrentFlowVelocity(local_idx,
                                                                  new_vel_edge);
    }
  });

  // remove empty particles
  removeEmptyParticles();
  makeParticlesConsistent();
}

void LiquidSimulator::transferLiquidFlowGrid() {
  const Scalar h_t = m_liquid_info.typical_flow_thickness;
  const Scalar dx = getDx();

  const int num_fluid = numParticles();

  if (num_fluid == 0) return;

  const int num_strand = m_strands.size();
  if (!num_strand) return;

  // prepare strands
  for_each(0, num_strand, [&](int strand_idx) {
    ElasticStrand* strand = m_strands[strand_idx];
    strand->dynamics().updateInterfaceSegments();
  });

  // prepare sorter
  m_particle_cells_single.sort(
      num_fluid, [&](int pidx, int& i, int& j, int& k) {
        Vec3x local_x = (m_x.segment<3>(pidx * 3) - m_grid_mincorner) / dx;
        i = (int)floor(local_x(0));
        j = (int)floor(local_x(1));
        k = (int)floor(local_x(2));
      });

  // for all grid node find inside vertices and register their extra flow volume
  // onto grid
  m_vertex_buckets.for_each_bucket_particles_colored([&](int pidx, int) {
    const int strand_idx = m_strands_global_local[pidx].first;
    const int local_idx = m_strands_global_local[pidx].second;

    // check if near to liquid
    const StrandDynamicTraits& dyn = m_strands[strand_idx]->dynamics();

    if (!dyn.isVertexNear(local_idx)) return;

    const Scalar len = m_strands[strand_idx]->getVoronoiLength(local_idx);
    const Scalar area = m_strands[strand_idx]->getCurrentFlowDOFArea(local_idx);

    const Scalar ra = m_strands[strand_idx]->getRadiusA(local_idx);
    const Scalar rb = m_strands[strand_idx]->getRadiusB(local_idx);
    // determain volume to be transferred
    Scalar vol;
    Scalar vol_on_flow;
    // check if inside the liquid
    if (dyn.isVertexInside(local_idx)) {
      // put all volume onto grid
      vol = area * len;
      vol_on_flow = 0.0;
    } else {
      // put extra volume onto grid
      const Scalar limit_area = M_PI * h_t * (h_t + ra + rb);

      vol = std::max(0.0, area - limit_area) * len;
      vol_on_flow = area * len - vol;
    }

    if (vol == 0.0) return;

    // find neighbor liquid particles
    const Vec3x pos = m_strands[strand_idx]->getVertex(local_idx);
    const Vec3i ipos = Vec3i((pos(0) - m_grid_mincorner(0)) / dx,
                             (pos(1) - m_grid_mincorner(1)) / dx,
                             (pos(2) - m_grid_mincorner(2)) / dx);

    const int cell_idx = m_particle_cells_single.bucket_index(ipos);

    const int num_verts = m_strands[strand_idx]->getNumVertices();

    const Vec3x u_s = dyn.getDisplacement(local_idx) / m_dt;
    Vec3x u_tau = Vec3x::Zero();
    int count = 0;
    if (local_idx > 0) {
      u_tau += m_strands[strand_idx]->getStepper()->getCurrentFlowVelocity()(
                   local_idx - 1) *
               m_strands[strand_idx]->getEdgeVector(local_idx - 1).normalized();
      count++;
    }
    if (local_idx < num_verts - 1) {
      u_tau += m_strands[strand_idx]->getStepper()->getCurrentFlowVelocity()(
                   local_idx) *
               m_strands[strand_idx]->getEdgeVector(local_idx).normalized();
      count++;
    }

    if (!count) return;
    u_tau /= (Scalar)count;

    Vec3x u_f = u_s + u_tau;

    const VecXx& c =
        m_strands[strand_idx]->getStepper()->flowComponents().segment(
            local_idx * m_liquid_info.num_components,
            m_liquid_info.num_components);

    Scalar vol_remain = vol;
    m_particle_cells_single.loop_neighbor_bucket_particles(
        cell_idx,
        [&](int npidx, int ncidx) {
          const Vec3x& np = m_x.segment<3>(npidx * 3);
          Vec3x diff = (pos - np) / dx;
          const Scalar w = N_kernel<2>(diff);
          Scalar vol_share = std::min(vol_remain, w * vol);

          Vec3x mom_share = u_f * vol_share;
          Vec3x mom_ori = m_v.segment<3>(npidx * 3) * m_rest_vol(npidx);
          VecXx color_share = c * vol_share;

          m_v.segment<3>(npidx * 3) =
              (mom_share + mom_ori) /
              std::max(1e-20, m_rest_vol(npidx) + vol_share);

          const VecXx color_ori =
              m_components.segment(npidx * m_liquid_info.num_components,
                                   m_liquid_info.num_components) *
              m_rest_vol(npidx);
          VecXx color_new = color_ori + color_share;
          make_gibbs_simplex(color_new);

          m_components.segment(npidx * m_liquid_info.num_components,
                               m_liquid_info.num_components) = color_new;

          m_vol(npidx) += vol_share;
          m_rest_vol(npidx) += vol_share;

          vol_remain -= vol_share;

          return vol_remain == 0.0;
        },
        2);

    vol_on_flow += vol_remain;

    // put back to flow
    m_strands[strand_idx]->setCurrentFlowDOFArea(local_idx, vol_on_flow / len);
  });

  makeParticlesComponentConsistent();
  makeParticlesConsistent();
}

void LiquidSimulator::applyGridBasedBodyFriction(const vector<VecXx>& lhs_x,
                                                 const vector<VecXx>& lhs_y,
                                                 const vector<VecXx>& lhs_z) {
  const Scalar dx = getDx();

  if (m_liquid_info.liquid_boundary_friction == 0.0) return;

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    const int num_nodes = getNumNodes(bucket_idx);

    auto project_vel = [&](vector<VecXx>& node_rhs_fluid,
                           const vector<VecXx>& node_solid_vel,
                           const vector<VecXx>& node_mass_fluid, const int ir) {
      for (int i = 0; i < num_nodes; ++i) {
        if (node_mass_fluid[bucket_idx][i] < 1e-20) continue;

        const Vec3x np = getNodePos(bucket_idx, i, ir);
        Vec3x n = Vec3x::Zero();
        const Scalar phi = interpolateValueAndGradient(
            n, np, m_node_solid_phi, m_grid_mincorner, 3.0 * dx);

        if (phi > dx * m_liquid_info.solid_shell_thickness) continue;

        n.normalize();

        Vec3x vel = get_future_velocity(np);
        Vec3x svel = get_solid_velocity(np);

        Vec3x rvel = vel - svel;
        Vec3x tvel = (Mat3x::Identity() - n * n.transpose()) * rvel;
        Vec3x nvel = n * n.transpose() * rvel;

        Vec3x proj_vel =
            clamp(1. - m_liquid_info.liquid_boundary_friction, 0.0, 1.0) *
                tvel +
            nvel + svel;

        // No Slip! Damn yeah no slip!
        node_rhs_fluid[bucket_idx][i] =
            proj_vel[ir] * node_mass_fluid[bucket_idx][i];
      }
    };

    project_vel(m_node_rhs_fluid_x, m_node_solid_vel_x, lhs_x, 0);
    project_vel(m_node_rhs_fluid_y, m_node_solid_vel_y, lhs_y, 1);
    project_vel(m_node_rhs_fluid_z, m_node_solid_vel_z, lhs_z, 2);
  });
}

void LiquidSimulator::conservativeResizeParticles(int num_particles) {
  m_x.conservativeResize(3 * num_particles);
  m_v.conservativeResize(3 * num_particles);
  m_m.conservativeResize(3 * num_particles);
  m_vol.conservativeResize(num_particles);
  m_rest_vol.conservativeResize(num_particles);
  m_radius.conservativeResize(num_particles);
  m_particle_group.resize(num_particles);
  m_J.conservativeResize(num_particles);

  m_B.conservativeResize(num_particles * 3, 3);
  m_Fe.conservativeResize(num_particles * 3, 3);
  m_Fe_plus.conservativeResize(num_particles * 3, 3);
  m_b.conservativeResize(num_particles * 3, 3);
  m_b_trial.conservativeResize(num_particles * 3, 3);
  m_b_plus.conservativeResize(num_particles * 3, 3);
  m_proj_func.conservativeResize(3 * num_particles);

  m_particle_nodes_x.resize(num_particles);
  m_particle_nodes_y.resize(num_particles);
  m_particle_nodes_z.resize(num_particles);
  m_particle_nodes_p.resize(num_particles);

  m_particle_weights.resize(num_particles);
  //        m_particle_grad_x.resize(num_particles);
  //        m_particle_grad_y.resize(num_particles);
  //        m_particle_grad_z.resize(num_particles);

  m_components.conservativeResize(num_particles * m_liquid_info.num_components);

  m_classifier.resize(num_particles);
  m_weakened.conservativeResize(num_particles);
}

int LiquidSimulator::getNumComponents() const {
  return m_liquid_info.num_components;
}

void LiquidSimulator::dripFromReservoir() {
  const Scalar dx = getDx();
  const int num_end_verts = m_strands.size() * 2;

  if (!num_end_verts) return;

  Vec3x bbx_min = Vec3x::Constant(std::numeric_limits<Scalar>::max());
  Vec3x bbx_max = Vec3x::Constant(-std::numeric_limits<Scalar>::max());
  for (int i = 0; i < num_end_verts; ++i) {
    int strand_idx = i / 2;
    int tip_idx = (m_strands[strand_idx]->getNumVertices() - 1) * (i & 1);

    Vec3x p = m_strands[strand_idx]->getVertex(tip_idx);

    bbx_min = Vec3x(std::min(p(0), bbx_min(0)), std::min(p(1), bbx_min(1)),
                    std::min(p(2), bbx_min(2)));
    bbx_max = Vec3x(std::max(p(0), bbx_max(0)), std::max(p(1), bbx_max(1)),
                    std::max(p(2), bbx_max(2)));
  }

  bbx_min = Vec3x(floor(bbx_min(0) / dx) * dx, floor(bbx_min(1) / dx) * dx,
                  floor(bbx_min(2) / dx) * dx);
  bbx_max = Vec3x(ceil(bbx_max(0) / dx) * dx, ceil(bbx_max(1) / dx) * dx,
                  ceil(bbx_max(2) / dx) * dx);

  bbx_min -= Vec3x::Constant(dx);
  bbx_max += Vec3x::Constant(dx);

  const Vec3x grid_size = bbx_max - bbx_min;

  const int ni = (int)ceil(grid_size(0) / dx);
  const int nj = (int)ceil(grid_size(1) / dx);
  const int nk = (int)ceil(grid_size(2) / dx);

  const int sni = (ni + m_num_nodes - 1) / m_num_nodes;
  const int snj = (nj + m_num_nodes - 1) / m_num_nodes;
  const int snk = (nk + m_num_nodes - 1) / m_num_nodes;

  struct Sample {
    Vec3x p;
    Vec3x v;
    VecXx c;
  };

  vector<vector<Sample> > bucket_samples(sni * snj * snk);

  Sorter sorter(ni, nj, nk);

  sorter.sort(num_end_verts, [&](int pidx, int& i, int& j, int& k) {
    int strand_idx = pidx / 2;
    int tip_idx = (m_strands[strand_idx]->getNumVertices() - 1) * (pidx & 1);

    Vec3x p = m_strands[strand_idx]->getVertex(tip_idx);

    i = (int)std::floor((p(0) - bbx_min(0)) / dx);
    j = (int)std::floor((p(1) - bbx_min(1)) / dx);
    k = (int)std::floor((p(2) - bbx_min(2)) / dx);
  });

  const int clampedInvCellMultiplier =
      (int)ceil(1.0 / m_liquid_info.particle_cell_multiplier);

  Sorter sorter_existed(ni * clampedInvCellMultiplier,
                        nj * clampedInvCellMultiplier,
                        nk * clampedInvCellMultiplier);
  const int num_particles = numParticles();
  const Scalar half_dx = dx / (Scalar)clampedInvCellMultiplier;
  sorter_existed.sort(num_particles, [&](int pidx, int& i, int& j, int& k) {
    i = (int)std::floor((m_x(pidx * 3 + 0) - bbx_min(0)) / half_dx);
    j = (int)std::floor((m_x(pidx * 3 + 1) - bbx_min(1)) / half_dx);
    k = (int)std::floor((m_x(pidx * 3 + 2) - bbx_min(2)) / half_dx);
  });

  const Scalar rad =
      defaultRadiusMultiplier() * dx * m_liquid_info.particle_cell_multiplier;
  const Scalar pvol = 4.0 / 3.0 * M_PI * rad * rad * rad;

  sorter.for_each_bucket_colored(
      [&](int cell_idx) {
        Scalar vol = 0.0;
        Vec3x vel = Vec3x::Zero();
        Vec3x pos = Vec3x::Zero();
        VecXx color = VecXx::Zero(m_liquid_info.num_components);

        // accumulate reservoir
        sorter.get_bucket(cell_idx, [&](int pidx) {
          int strand_idx = pidx / 2;
          int tip_idx =
              (m_strands[strand_idx]->getNumVertices() - 1) * (pidx & 1);
          int tip_edge_idx =
              (m_strands[strand_idx]->getNumEdges() - 1) * (pidx & 1);

          Vec2x& reservoir = m_strands[strand_idx]->getReservoir();
          Scalar w = reservoir[pidx & 1];
          vol += w;

          const Scalar ut =
              m_strands[strand_idx]->getStepper()->getCurrentFlowVelocity()(
                  tip_edge_idx);
          const Vec3x ev =
              m_strands[strand_idx]->getEdgeVector(tip_edge_idx).normalized();

          vel +=
              (m_strands[strand_idx]->dynamics().getDisplacements().segment<3>(
                   tip_idx * 4) /
                   m_dt +
               ut * ev) *
              w;
          pos += m_strands[strand_idx]->getVertex(tip_idx) * w;
          color +=
              m_strands[strand_idx]->getStepper()->flowComponents().segment(
                  tip_idx * m_liquid_info.num_components,
                  m_liquid_info.num_components) *
              w;
        });

        if (vol < pvol) return;

        int num_release = std::min((int)floor(vol / pvol),
                                   (int)sphere_pattern::max_vector_length);
        Scalar actual_released_vol = pvol * num_release;

        vel /= vol;
        pos /= vol;
        Scalar cs = color.sum();
        if (cs == 0.0) {
          color.setZero();
          color(0) = 1.0;
        } else {
          color /= cs;
        }

        Vec3i cell_handle = sorter.bucket_handle(cell_idx);
        Vec3i bucket_handle =
            Vec3i(cell_handle(0) / m_num_nodes, cell_handle(1) / m_num_nodes,
                  cell_handle(2) / m_num_nodes);
        const int bucket_idx = bucket_handle(2) * sni * snj +
                               bucket_handle(1) * sni + bucket_handle(0);

        const Scalar radius_released_vol =
            pow(actual_released_vol * 0.75 / M_PI, 1.0 / 3.0);

        // generate particles
        Mat3x M = Mat3x::Random();
        Mat3x Q, R;
        QRDecompose<Scalar, 3>(M, Q, R);

        for (int i = 0; i < num_release; ++i) {
          Sample s;
          Vec3x actual_pos =
              pos + Q * m_sphere_pattern[num_release].segment<3>(i * 3) *
                        radius_released_vol;

          // check if pos available (not blocked by existing particles)
          Vec3i ipos_existed =
              Vec3i((int)std::floor((actual_pos(0) - bbx_min(0)) / half_dx),
                    (int)std::floor((actual_pos(1) - bbx_min(1)) / half_dx),
                    (int)std::floor((actual_pos(2) - bbx_min(2)) / half_dx));

          if (sorter_existed.has_bucket(ipos_existed)) {
            const int cell_idx_existed =
                sorter_existed.bucket_index(ipos_existed);

            // remove the particle if occupied by existing particles
            if (sorter_existed.get_bucket_size(cell_idx_existed) > 0) {
              actual_released_vol -= pvol;
              continue;
            }
          }

          // pos is available
          s.v = vel;
          s.p = actual_pos;
          s.c = color;
          bucket_samples[bucket_idx].push_back(s);
        }

        // remove from reservoir
        const Scalar prop = 1.0 - actual_released_vol / vol;
        sorter.get_bucket(cell_idx, [&](int pidx) {
          int strand_idx = pidx / 2;
          int tip_idx =
              (m_strands[strand_idx]->getNumVertices() - 1) * (pidx & 1);
          int tip_edge_idx =
              (m_strands[strand_idx]->getNumEdges() - 1) * (pidx & 1);

          Vec2x& reservoir = m_strands[strand_idx]->getReservoir();
          reservoir[pidx & 1] *= prop;
        });
      },
      m_num_nodes);

  vector<Sample> new_samples;
  parallel_concatenate(new_samples, bucket_samples);

  const int num_new_samples = new_samples.size();
  if (!num_new_samples) return;

  {
    LockGuard lock(m_particle_mutex);

    const int num_old_parts = numParticles();
    conservativeResizeParticles(num_old_parts + num_new_samples);

    for_each(0, num_new_samples, [&](int i) {
      const int part_idx = num_old_parts + i;
      m_x.segment<3>(part_idx * 3) = new_samples[i].p;
      m_v.segment<3>(part_idx * 3) = new_samples[i].v;
      m_m.segment<3>(part_idx * 3)
          .setConstant(pvol * getDensity(new_samples[i].c));
      m_vol(part_idx) = m_rest_vol(part_idx) = pvol;
      m_radius(part_idx) = rad;
      m_particle_group[part_idx] = 0;
      m_B.block<3, 3>(part_idx * 3, 0).setZero();
      m_Fe.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_b.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_b_trial.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_b_plus.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_Fe_plus.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_proj_func.segment<3>(part_idx * 3) = Vec3x(1, 0, 0);
      m_J(part_idx) = 1.0;
      m_weakened[part_idx] = false;
      m_components.segment(part_idx * m_liquid_info.num_components,
                           m_liquid_info.num_components) = new_samples[i].c;
    });
  }
}

bool LiquidSimulator::useConstantDrag() const {
  return m_liquid_info.use_constant_drag;
}

Scalar LiquidSimulator::getDragInsulation() const {
  return m_liquid_info.geometric_drag_insulation;
}

void LiquidSimulator::sampleLiquidDistanceFields(Scalar cur_time) {
  LockGuard lock(m_particle_mutex);

  int num_group = (int)m_sources.size();

  const Scalar dx = getDx();  // we use denser dx to prevent penetration

  for (int igroup = 0; igroup < num_group; ++igroup) {
    Vec3x shooting_vel = Vec3x::Zero();

    if (m_sources[igroup].usage != DFU_SOURCE ||
        !m_sources[igroup].check_durations(
            cur_time, m_shooting_vol_accum[igroup], shooting_vel))
      continue;

    VecXx additional_pos;

    m_sources[igroup].resample_internal(
        m_fields, dx * m_liquid_info.particle_cell_multiplier, m_x,
        additional_pos,
        (int)ceil((Scalar)m_num_nodes /
                  m_liquid_info.particle_cell_multiplier));

    const int df_index = numParticles();
    const int df_size = additional_pos.size() / 3;

    if (df_size == 0) continue;

    conservativeResizeParticles(df_index + df_size);

    const Scalar rad =
        defaultRadiusMultiplier() * dx * m_liquid_info.particle_cell_multiplier;

    const Scalar pvol = 4.0 / 3.0 * M_PI * rad * rad * rad;

    m_shooting_vol_accum[igroup] += pvol * (Scalar)df_size;

    for_each(0, df_size, [&](int i) {
      const int part_idx = df_index + i;
      Vec3x pos = additional_pos.segment<3>(i * 3);
      Vec3x vel = Vec3x::Zero();
      m_sources[igroup].compute_phi_vel(pos, vel);

      m_x.segment<3>(part_idx * 3) = pos;
      m_v.segment<3>(part_idx * 3) = vel + shooting_vel;

      VecXx color = VecXx::Zero(m_liquid_info.num_components);
      color(m_sources[igroup].color_index) = 1.0;

      m_m.segment<3>(part_idx * 3).setConstant(pvol * getDensity(color));
      m_vol(part_idx) = m_rest_vol(part_idx) = pvol;
      m_radius(part_idx) = rad;
      m_particle_group[part_idx] = igroup;
      m_B.block<3, 3>(part_idx * 3, 0).setZero();
      m_Fe.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_b.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_b_trial.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_b_plus.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_Fe_plus.block<3, 3>(part_idx * 3, 0).setIdentity();
      m_proj_func.segment<3>(part_idx * 3) = Vec3x(1, 0, 0);
      m_J(part_idx) = 1.0;
      m_weakened[part_idx] = false;
      m_components.segment(part_idx * m_liquid_info.num_components,
                           m_liquid_info.num_components) = color;
    });
  }

  // check_isnan("x_after_sampled", m_x);
  // check_isnan("v_after_sampled", m_v);
  // check_isnan("m_after_sampled", m_m);
}

void LiquidSimulator::updateParticleBoundingBox() {
  const int numParts = numParticles();

  if (numParts == 0 && !m_strands.size()) {
    m_bbx_min = m_bbx_max = Vec3x::Zero();
  } else {
    Vec3x bbmin = Vec3x::Constant(1e+20);
    Vec3x bbmax = Vec3x::Constant(-1e+20);

    bbmin = reduction((Vec3x*)m_x.data(), numParticles(), bbmin,
                      [](Vec3x x, Vec3x y) -> Vec3x {
                        return Vec3x(std::min(x(0), y(0)), std::min(x(1), y(1)),
                                     std::min(x(2), y(2)));
                      });

    bbmax = reduction((Vec3x*)m_x.data(), numParticles(), bbmax,
                      [](Vec3x x, Vec3x y) -> Vec3x {
                        return Vec3x(std::max(x(0), y(0)), std::max(x(1), y(1)),
                                     std::max(x(2), y(2)));
                      });

    const Scalar dx = m_bucket_size * 2.0;

    std::vector<Vec3x> bbmins(m_strands.size(), bbmin);
    std::vector<Vec3x> bbmaxs(m_strands.size(), bbmax);

    const int num_strands = m_strands.size();

    for_each(0, num_strands, [&](int k) {
      const ElasticStrand* es = m_strands[k];
      const int num_verts = es->getNumVertices();
      for (int i = 0; i < num_verts; ++i) {
        Vec3x p = es->getVertex(i);
        bbmins[k](0) = std::min(bbmins[k](0), p(0));
        bbmins[k](1) = std::min(bbmins[k](1), p(1));
        bbmins[k](2) = std::min(bbmins[k](2), p(2));
        bbmaxs[k](0) = std::max(bbmaxs[k](0), p(0));
        bbmaxs[k](1) = std::max(bbmaxs[k](1), p(1));
        bbmaxs[k](2) = std::max(bbmaxs[k](2), p(2));
      }
    });

    for (int i = 0; i < num_strands; ++i) {
      bbmin(0) = std::min(bbmin(0), bbmins[i](0));
      bbmin(1) = std::min(bbmin(1), bbmins[i](1));
      bbmin(2) = std::min(bbmin(2), bbmins[i](2));
      bbmax(0) = std::max(bbmax(0), bbmaxs[i](0));
      bbmax(1) = std::max(bbmax(1), bbmaxs[i](1));
      bbmax(2) = std::max(bbmax(2), bbmaxs[i](2));
    }

    m_bbx_min = Vec3x(floor(bbmin(0) / dx) * dx, floor(bbmin(1) / dx) * dx,
                      floor(bbmin(2) / dx) * dx);
    m_bbx_max = Vec3x(ceil(bbmax(0) / dx) * dx, ceil(bbmax(1) / dx) * dx,
                      ceil(bbmax(2) / dx) * dx);
  }

  std::cout << "[bounding box min: " << m_bbx_min(0) << ", " << m_bbx_min(1)
            << ", " << m_bbx_min(2) << "; max: " << m_bbx_max(0) << ", "
            << m_bbx_max(1) << ", " << m_bbx_max(2) << "]" << std::endl;
}

void LiquidSimulator::rebucketizeParticles() {
  Scalar dx = getDx();

  const Scalar extra_border = 3.0;

  Vec3x content_size = m_bbx_max - m_bbx_min +
                       Vec3x::Constant(m_bucket_size * extra_border * 2.0);

  Vec3i grid_num_cells =
      Vec3i(std::max(1, (int)std::ceil(content_size(0) / dx)),
            std::max(1, (int)std::ceil(content_size(1) / dx)),
            std::max(1, (int)std::ceil(content_size(2) / dx)));

  Vec3x grid_size =
      Vec3x((Scalar)grid_num_cells[0] * dx, (Scalar)grid_num_cells[1] * dx,
            (Scalar)grid_num_cells[2] * dx);

  m_grid_mincorner = m_bbx_min - Vec3x::Constant(m_bucket_size * extra_border);

  Vec3i num_buckets =
      Vec3i(std::max(1, (int)std::ceil(grid_size(0) / m_bucket_size)),
            std::max(1, (int)std::ceil(grid_size(1) / m_bucket_size)),
            std::max(1, (int)std::ceil(grid_size(2) / m_bucket_size)));

  m_particle_buckets.resize(num_buckets(0), num_buckets(1), num_buckets(2));
  m_particle_cells_single.resize(grid_num_cells(0), grid_num_cells(1),
                                 grid_num_cells(2));

  m_particle_buckets.sort(
      numParticles(), [&](int pidx, int& i, int& j, int& k) {
        i = (int)std::floor((m_x(pidx * 3 + 0) - m_grid_mincorner(0)) /
                            m_bucket_size);
        j = (int)std::floor((m_x(pidx * 3 + 1) - m_grid_mincorner(1)) /
                            m_bucket_size);
        k = (int)std::floor((m_x(pidx * 3 + 2) - m_grid_mincorner(2)) /
                            m_bucket_size);
      });

  const int total_buckets = m_particle_buckets.size();

  m_bucket_activated.assign(total_buckets, 0U);

  m_vertex_buckets.resize(num_buckets(0), num_buckets(1), num_buckets(2));
  m_vertex_buckets.sort(
      m_strands_global_local.size(), [&](int pidx, int& i, int& j, int& k) {
        const pair<int, int>& local_idx = m_strands_global_local[pidx];
        const ElasticStrand* strand = m_strands[local_idx.first];
        const Vec3x pos = strand->getVertex(local_idx.second);
        i = (int)std::floor((pos(0) - m_grid_mincorner(0)) / m_bucket_size);
        j = (int)std::floor((pos(1) - m_grid_mincorner(1)) / m_bucket_size);
        k = (int)std::floor((pos(2) - m_grid_mincorner(2)) / m_bucket_size);
      });

  m_mesh_vertex_buckets.resize(num_buckets(0), num_buckets(1), num_buckets(2));

  std::cout << "[num buckets: " << num_buckets(0) << ", " << num_buckets(1)
            << ", " << num_buckets(2) << "]" << std::endl;
}

void LiquidSimulator::preAllocateNodes() {
  // mark buckets as active
  const int num_buckets = m_particle_buckets.size();

  m_node_pos.resize(num_buckets);

  m_node_pressure_neighbors.resize(num_buckets);
  m_node_pp_neighbors.resize(num_buckets);
  m_node_index_pressure_x.resize(num_buckets);
  m_node_index_pressure_y.resize(num_buckets);
  m_node_index_pressure_z.resize(num_buckets);

  m_node_index_solid_phi_x.resize(num_buckets);
  m_node_index_solid_phi_y.resize(num_buckets);
  m_node_index_solid_phi_z.resize(num_buckets);
}

int LiquidSimulator::getDefaultNumNodes() const { return m_num_nodes; }

Scalar LiquidSimulator::get_dense_liquid_signed_distance(
    const Vec3x& pos) const {
  return interpolateValue(pos, m_node_signed_distance, m_grid_mincorner, 0.0,
                          m_liquid_info.signed_distance_multiplier);
}

Scalar LiquidSimulator::getCellSize() const {
  return m_bucket_size / (Scalar)m_num_nodes;
}

const vector<VecXx>& LiquidSimulator::getNodeSignedDistance() const {
  return m_node_signed_distance;
}

bool LiquidSimulator::computeSignedDistance() {
  const int num_buckets = m_particle_buckets.size();
  if (!num_buckets) return true;

  LockGuard lock(m_grid_mutex);

  const Scalar dx = getDx();
  const Scalar dr = dx / (Scalar)m_liquid_info.signed_distance_multiplier;

  const int dni = m_particle_buckets.ni * m_num_nodes *
                  m_liquid_info.signed_distance_multiplier;
  const int dnj = m_particle_buckets.nj * m_num_nodes *
                  m_liquid_info.signed_distance_multiplier;
  const int dnk = m_particle_buckets.nk * m_num_nodes *
                  m_liquid_info.signed_distance_multiplier;

  //       std::cout << "[computeSignedDistance 0]" << std::endl;

  std::cout << "[csd num nodes: " << dni << ", " << dnj << ", " << dnk << "]"
            << std::endl;

  m_particle_cells.resize(dni, dnj, dnk);

  m_particle_cells.sort(numParticles(), [&](int pidx, int& i, int& j, int& k) {
    i = (int)std::floor((m_x(pidx * 3 + 0) - m_grid_mincorner(0)) / dr);
    j = (int)std::floor((m_x(pidx * 3 + 1) - m_grid_mincorner(1)) / dr);
    k = (int)std::floor((m_x(pidx * 3 + 2) - m_grid_mincorner(2)) / dr);
  });

  allocateNodeVectors(m_node_signed_distance,
                      m_liquid_info.signed_distance_multiplier *
                          m_liquid_info.signed_distance_multiplier *
                          m_liquid_info.signed_distance_multiplier,
                      std::numeric_limits<Scalar>::max());

  //       std::cout << "[computeSignedDistance 1]" << std::endl;

  const int num_nodes_dist =
      m_num_nodes * m_liquid_info.signed_distance_multiplier;

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    const Vec3x base_np = getNodePos(bucket_idx, 0, 4);
    const Vec3i base_handle =
        m_particle_buckets.bucket_handle(bucket_idx) * num_nodes_dist;

    VecXx& bucket_signed_distance = m_node_signed_distance[bucket_idx];

    for (int k = 0; k < num_nodes_dist; ++k)
      for (int j = 0; j < num_nodes_dist; ++j)
        for (int i = 0; i < num_nodes_dist; ++i) {
          const int local_node_idx =
              k * num_nodes_dist * num_nodes_dist + j * num_nodes_dist + i;

          const Vec3x np = base_np + Vec3x(i, j, k) * dr;
          const Vec3i node_handle = base_handle + Vec3i(i, j, k);
          const int node_idx = m_particle_cells.bucket_index(node_handle);

          // init with solid SDF here
          bucket_signed_distance[local_node_idx] =
              DistanceFieldObject::computePhi(np, m_fields);

          m_particle_cells.loop_neighbor_bucket_particles(
              node_idx,
              [&](int pidx, int) -> bool {
                const Scalar dist = (m_x.segment<3>(pidx * 3) - np).norm();
                if (dist >= dx) return false;

                if (dist < bucket_signed_distance[local_node_idx]) {
                  bucket_signed_distance[local_node_idx] = dist;
                }
                return false;
              },
              2, 1, 2, 1, 2, 1);

          bucket_signed_distance[local_node_idx] -= dr;
        }
  });

  //       std::cout << "[computeSignedDistance 2]" << std::endl;

  vector<vector<Vec3x> > interface_points(num_buckets);

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    vector<Vec3x>& bucket_ip = interface_points[bucket_idx];
    bucket_ip.reserve(num_nodes_dist * num_nodes_dist * 3);

    Vec3i bucket_handle = m_particle_buckets.bucket_handle(bucket_idx);
    const Vec3x base_np = getNodePos(bucket_idx, 0, 4);

    for (int k = 0; k < num_nodes_dist; ++k)
      for (int j = 0; j < num_nodes_dist; ++j)
        for (int i = 0; i < num_nodes_dist; ++i) {
          Vec3i cur_node_handle = Vec3i(i, j, k);
          const Vec3x cur_np = base_np + Vec3x(i, j, k) * dr;

          const int cur_node_local_index =
              k * num_nodes_dist * num_nodes_dist + j * num_nodes_dist + i;

          const Scalar& cur_node_dist =
              m_node_signed_distance[bucket_idx][cur_node_local_index];

          for (int r = 0; r < 3; ++r) {
            Vec3i next_node_handle = cur_node_handle + Vec3i::Unit(r);
            Vec3i next_bucket_handle = bucket_handle;
            if (next_node_handle(r) >= num_nodes_dist) {
              next_bucket_handle(r)++;
              next_node_handle(r) -= num_nodes_dist;
            }

            if (next_bucket_handle(r) >= m_particle_buckets.dim_size(r))
              continue;

            const int next_bucket_idx =
                m_particle_buckets.bucket_index(next_bucket_handle);
            if (!m_bucket_activated[next_bucket_idx]) continue;

            const int next_node_local_index =
                next_node_handle(2) * num_nodes_dist * num_nodes_dist +
                next_node_handle(1) * num_nodes_dist + next_node_handle(0);

            const Scalar& next_node_dist =
                m_node_signed_distance[next_bucket_idx][next_node_local_index];

            if (next_node_dist * cur_node_dist <= 0.0) {
              Vec3x spos = cur_np - cur_node_dist * Vec3x::Unit(r);
              bucket_ip.push_back(spos);
            }
          }
        }
  });

  //        std::cout << "[computeSignedDistance 3]" << std::endl;

  vector<int> bucket_num_ips(num_buckets);
  // parallel concatenate
  for_each(0, num_buckets, [&](int bucket_idx) {
    bucket_num_ips[bucket_idx] = interface_points[bucket_idx].size();
  });

  std::partial_sum(bucket_num_ips.begin(), bucket_num_ips.end(),
                   bucket_num_ips.begin());
  const int num_ips = bucket_num_ips[num_buckets - 1];

  m_interface_points.resize(num_ips * 3);
  if (num_ips == 0) return true;

  //       std::cout << "[computeSignedDistance 4]" << std::endl;

  for_each(0, num_buckets, [&](int bucket_idx) {
    int base_idx = (bucket_idx == 0) ? 0 : bucket_num_ips[bucket_idx - 1];
    const int num_ips_local = interface_points[bucket_idx].size();
    for (int i = 0; i < num_ips_local; ++i) {
      m_interface_points.segment<3>((base_idx + i) * 3) =
          interface_points[bucket_idx][i];
    }
  });

  //       std::cout << "[computeSignedDistance 5]" << std::endl;

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    const int num_cells = m_node_signed_distance[bucket_idx].size();
    for (int i = 0; i < num_cells; ++i) {
      m_node_signed_distance[bucket_idx][i] =
          sgn(m_node_signed_distance[bucket_idx][i]) *
          std::numeric_limits<Scalar>::max();
    }
  });

  //       std::cout << "[computeSignedDistance 6]" << std::endl;

  m_interface_buckets.resize(m_particle_buckets.ni, m_particle_buckets.nj,
                             m_particle_buckets.nk);
  m_interface_buckets.sort(num_ips, [&](int pidx, int& i, int& j, int& k) {
    i = (int)std::floor(
        (m_interface_points(pidx * 3 + 0) - m_grid_mincorner(0)) /
        m_bucket_size);
    j = (int)std::floor(
        (m_interface_points(pidx * 3 + 1) - m_grid_mincorner(1)) /
        m_bucket_size);
    k = (int)std::floor(
        (m_interface_points(pidx * 3 + 2) - m_grid_mincorner(2)) /
        m_bucket_size);
  });

  //       std::cout << "[computeSignedDistance 7]" << std::endl;

  m_interface_buckets.for_each_bucket_particles_colored_randomized(
      [&](int pidx, int bucket_idx) {
        const Vec3x& ip_pos = m_interface_points.segment<3>(pidx * 3);
        const Vec3i global_cell_handle =
            Vec3i((int)std::floor((ip_pos(0) - m_grid_mincorner(0)) / dr),
                  (int)std::floor((ip_pos(1) - m_grid_mincorner(1)) / dr),
                  (int)std::floor((ip_pos(2) - m_grid_mincorner(2)) / dr));

        const Vec3i bucket_handle =
            m_interface_buckets.bucket_handle(bucket_idx);
        const Vec3i local_cell_handle =
            Vec3i(global_cell_handle(0) - bucket_handle(0) * num_nodes_dist,
                  global_cell_handle(1) - bucket_handle(1) * num_nodes_dist,
                  global_cell_handle(2) - bucket_handle(2) * num_nodes_dist);

        const Vec3x np_base =
            getNodePos(bucket_idx, 0, 4) + Vec3x(local_cell_handle(0) * dr,
                                                 local_cell_handle(1) * dr,
                                                 local_cell_handle(2) * dr);

        for (int k = -4; k <= 4; ++k)
          for (int j = -4; j <= 4; ++j)
            for (int i = -4; i <= 4; ++i) {
              const Vec3x np = np_base + Vec3x(i, j, k) * dr;
              const Scalar dist = (np - ip_pos).norm();
              if (dist > 4.0 * dr) continue;

              Vec3i neighbor_cell_handle = local_cell_handle + Vec3i(i, j, k);
              Vec3i neighbor_bucket_handle = bucket_handle;
              for (int r = 0; r < 3; ++r) {
                if (neighbor_cell_handle(r) < 0) {
                  neighbor_bucket_handle(r)--;
                  neighbor_cell_handle(r) += num_nodes_dist;
                } else if (neighbor_cell_handle(r) >= num_nodes_dist) {
                  neighbor_bucket_handle(r)++;
                  neighbor_cell_handle(r) -= num_nodes_dist;
                }
              }

              if (!m_interface_buckets.has_bucket(neighbor_bucket_handle))
                continue;
              const int neigh_bucket_idx =
                  m_interface_buckets.bucket_index(neighbor_bucket_handle);
              if (!m_bucket_activated[neigh_bucket_idx]) continue;

              const int neigh_node_idx =
                  neighbor_cell_handle(2) * num_nodes_dist * num_nodes_dist +
                  neighbor_cell_handle(1) * num_nodes_dist +
                  neighbor_cell_handle(0);

              if (dist < fabs(m_node_signed_distance[neigh_bucket_idx]
                                                    [neigh_node_idx])) {
                m_node_signed_distance[neigh_bucket_idx][neigh_node_idx] =
                    sgn(m_node_signed_distance[neigh_bucket_idx]
                                              [neigh_node_idx]) *
                    dist;
              }
            }
      },
      3);

  //       std::cout << "[computeSignedDistance 8]" << std::endl;

  return true;
}

VecXx LiquidSimulator::getComponents(int i) const {
  return m_components.segment(i * m_liquid_info.num_components,
                              m_liquid_info.num_components);
}

void LiquidSimulator::findNodes(const Sorter& buckets,
                                std::vector<Mat27x2i>& particle_nodes,
                                const Vec3x& offset,
                                const std::function<Vec3x(int)> get_pos,
                                bool activate) {
  const Scalar dx = getDx();

  for (Mat27x2i& term : particle_nodes) {
    term.setConstant(-1);
  }

  // make connection for particles
  buckets.for_each_bucket_colored([&](int bucket_idx) {
    Vec3i bucket_handle = buckets.bucket_handle(bucket_idx);

    Vec3x cell_local_corner = Vec3x((Scalar)bucket_handle(0) * m_bucket_size,
                                    (Scalar)bucket_handle(1) * m_bucket_size,
                                    (Scalar)bucket_handle(2) * m_bucket_size) +
                              m_grid_mincorner + offset * dx;

    buckets.get_bucket(bucket_idx, [&](int pidx) {
      Mat27x2i& indices = particle_nodes[pidx];

      Vec3x local_pos = (get_pos(pidx) - cell_local_corner) / dx;
      Vec3i ilocal_pos =
          Vec3i((int)floor(local_pos(0)), (int)floor(local_pos(1)),
                (int)floor(local_pos(2)));

      Vec3x local_frac = frac<Scalar, 3, 1>(local_pos);

      const int klow = local_frac(2) > 0.5 ? 0 : -1;
      const int jlow = local_frac(1) > 0.5 ? 0 : -1;
      const int ilow = local_frac(0) > 0.5 ? 0 : -1;
      const int khigh = klow + 2;
      const int jhigh = jlow + 2;
      const int ihigh = ilow + 2;

      for (int k = klow; k <= khigh; ++k)
        for (int j = jlow; j <= jhigh; ++j)
          for (int i = ilow; i <= ihigh; ++i) {
            Vec3i cell_local_idx = ilocal_pos + Vec3i(i, j, k);
            Vec3i node_bucket_handle = bucket_handle;

            bool out_of_range = false;
            for (int r = 0; r < 3; ++r) {
              while (cell_local_idx(r) < 0) {
                node_bucket_handle(r)--;
                cell_local_idx(r) += m_num_nodes;
              }
              while (cell_local_idx(r) >= m_num_nodes) {
                node_bucket_handle(r)++;
                cell_local_idx(r) -= m_num_nodes;
              }

              assert(cell_local_idx(r) >= 0 && cell_local_idx(r) < m_num_nodes);
              if (node_bucket_handle(r) < 0 ||
                  node_bucket_handle(r) >= buckets.dim_size(r)) {
                out_of_range = true;
                break;
              }
            }

            if (out_of_range) {
              int nidx = (k - klow) * (3 * 3) + (j - jlow) * 3 + (i - ilow);
              indices(nidx, 0) = indices(nidx, 1) = -1;
            } else {
              int node_bucket_idx = buckets.bucket_index(node_bucket_handle);
              if (activate) {
                m_bucket_activated[node_bucket_idx] = 1U;
              }

              int cell_idx = cell_local_idx(2) * m_num_nodes * m_num_nodes +
                             cell_local_idx(1) * m_num_nodes +
                             cell_local_idx(0);

              int nidx = (k - klow) * (3 * 3) + (j - jlow) * 3 + (i - ilow);
              indices(nidx, 0) = node_bucket_idx;
              indices(nidx, 1) = cell_idx;
            }
          }
    });
  });
}

void LiquidSimulator::expandFluidNodesMarked(int layers) {
  auto check_bucket = [&](const Vec3i& bucket_handle,
                          const std::vector<unsigned char>& activated) {
    for (int t = -1; t <= 1; ++t)
      for (int s = -1; s <= 1; ++s)
        for (int r = -1; r <= 1; ++r) {
          if (t == 0 && s == 0 && r == 0) continue;

          Vec3i cur_bucket_handle = bucket_handle + Vec3i(r, s, t);
          if (cur_bucket_handle(0) < 0 ||
              cur_bucket_handle(0) >= m_particle_buckets.ni ||
              cur_bucket_handle(1) < 0 ||
              cur_bucket_handle(1) >= m_particle_buckets.nj ||
              cur_bucket_handle(2) < 0 ||
              cur_bucket_handle(2) >= m_particle_buckets.nk) {
            continue;
          }

          const int nbidx = m_particle_buckets.bucket_index(cur_bucket_handle);
          if (activated[nbidx]) {
            return true;
          }
        }

    return false;
  };

  for (int iLayer = 0; iLayer < layers; ++iLayer) {
    std::vector<unsigned char> activated_backup = m_bucket_activated;

    m_particle_buckets.for_each_bucket([&](int bucket_idx) {
      const Vec3i bucket_handle = m_particle_buckets.bucket_handle(bucket_idx);

      // ignore bucket already activated or no neighbor bucket activated
      if (activated_backup[bucket_idx] ||
          !check_bucket(bucket_handle, activated_backup))
        return;

      // Activate bucket that has activated neighbor buckets.
      m_bucket_activated[bucket_idx] = true;
    });
  }
}

void LiquidSimulator::generateNodes() {
  const Scalar dx = getDx();

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    // ignore inactivated buckets
    if (!m_bucket_activated[bucket_idx]) return;

    Vec3i bucket_handle = m_particle_buckets.bucket_handle(bucket_idx);

    Vec3x cell_local_corner = Vec3x((Scalar)bucket_handle(0) * m_bucket_size,
                                    (Scalar)bucket_handle(1) * m_bucket_size,
                                    (Scalar)bucket_handle(2) * m_bucket_size) +
                              m_grid_mincorner;

    // otherwise generate nodes to fill the bucket
    const int count = m_num_nodes * m_num_nodes * m_num_nodes;
    if ((int)m_node_pos[bucket_idx].size() != count * 3)
      m_node_pos[bucket_idx].resize(count * 3);

    VecXx& bucket_node_pos = m_node_pos[bucket_idx];

    for (int k = 0; k < m_num_nodes; ++k)
      for (int j = 0; j < m_num_nodes; ++j)
        for (int i = 0; i < m_num_nodes; ++i) {
          int node_idx = k * m_num_nodes * m_num_nodes + j * m_num_nodes + i;
          bucket_node_pos.segment<3>(node_idx * 3) =
              cell_local_corner + Vec3x(i, j, k) * dx;
        }

    VecXi& bucket_node_idxp_x = m_node_index_pressure_x[bucket_idx];
    VecXi& bucket_node_idx_solid_phi_x = m_node_index_solid_phi_x[bucket_idx];
    bucket_node_idxp_x.resize(count * 4);
    bucket_node_idxp_x.setConstant(-1);
    bucket_node_idx_solid_phi_x.resize(count * 8);
    bucket_node_idx_solid_phi_x.setConstant(-1);

    VecXi& bucket_node_idxp_y = m_node_index_pressure_y[bucket_idx];
    VecXi& bucket_node_idx_solid_phi_y = m_node_index_solid_phi_y[bucket_idx];
    bucket_node_idxp_y.resize(count * 4);
    bucket_node_idxp_y.setConstant(-1);
    bucket_node_idx_solid_phi_y.resize(count * 8);
    bucket_node_idx_solid_phi_y.setConstant(-1);

    VecXi& bucket_node_idxp_z = m_node_index_pressure_z[bucket_idx];
    VecXi& bucket_node_idx_solid_phi_z = m_node_index_solid_phi_z[bucket_idx];
    bucket_node_idxp_z.resize(count * 4);
    bucket_node_idxp_z.setConstant(-1);
    bucket_node_idx_solid_phi_z.resize(count * 8);
    bucket_node_idx_solid_phi_z.setConstant(-1);
  });
}

void LiquidSimulator::connectSolidPhiNodes() {
  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    Vec3i bucket_handle = m_particle_buckets.bucket_handle(bucket_idx);

    VecXi& bucket_node_idx_solid_phi_x = m_node_index_solid_phi_x[bucket_idx];
    VecXi& bucket_node_idx_solid_phi_y = m_node_index_solid_phi_y[bucket_idx];
    VecXi& bucket_node_idx_solid_phi_z = m_node_index_solid_phi_z[bucket_idx];

    if (!m_bucket_activated[bucket_idx]) return;

    for (int k = 0; k < m_num_nodes; ++k)
      for (int j = 0; j < m_num_nodes; ++j)
        for (int i = 0; i < m_num_nodes; ++i) {
          int node_idx = k * m_num_nodes * m_num_nodes + j * m_num_nodes + i;

          for (int r = 0; r < 2; ++r)
            for (int s = 0; s < 2; ++s) {
              Vec3i node_bucket_handle_x = bucket_handle;

              Vec3i sphi_local_idx = Vec3i(i, j + s, k + r);

              if (sphi_local_idx(1) >= m_num_nodes) {
                node_bucket_handle_x(1)++;
                sphi_local_idx(1) -= m_num_nodes;
              }

              if (sphi_local_idx(2) >= m_num_nodes) {
                node_bucket_handle_x(2)++;
                sphi_local_idx(2) -= m_num_nodes;
              }

              if (node_bucket_handle_x(0) < 0 ||
                  node_bucket_handle_x(0) >= m_particle_buckets.dim_size(0) ||
                  node_bucket_handle_x(1) < 0 ||
                  node_bucket_handle_x(1) >= m_particle_buckets.dim_size(1) ||
                  node_bucket_handle_x(2) < 0 ||
                  node_bucket_handle_x(2) >= m_particle_buckets.dim_size(2)) {
                bucket_node_idx_solid_phi_x(node_idx * 8 + (r * 2 + s) * 2 +
                                            0) = -1;
                bucket_node_idx_solid_phi_x(node_idx * 8 + (r * 2 + s) * 2 +
                                            1) = -1;
              } else {
                const int sphi_idx =
                    sphi_local_idx(2) * m_num_nodes * m_num_nodes +
                    sphi_local_idx(1) * m_num_nodes + sphi_local_idx(0);

                int nb_bucket_idx =
                    m_particle_buckets.bucket_index(node_bucket_handle_x);
                if (!m_bucket_activated[nb_bucket_idx]) {
                  bucket_node_idx_solid_phi_x(node_idx * 8 + (r * 2 + s) * 2 +
                                              0) = -1;
                  bucket_node_idx_solid_phi_x(node_idx * 8 + (r * 2 + s) * 2 +
                                              1) = -1;
                } else {
                  bucket_node_idx_solid_phi_x(node_idx * 8 + (r * 2 + s) * 2 +
                                              0) = nb_bucket_idx;
                  bucket_node_idx_solid_phi_x(node_idx * 8 + (r * 2 + s) * 2 +
                                              1) = sphi_idx;
                }
              }
            }

          for (int r = 0; r < 2; ++r)
            for (int s = 0; s < 2; ++s) {
              Vec3i node_bucket_handle_y = bucket_handle;

              Vec3i sphi_local_idx = Vec3i(i + r, j, k + s);

              if (sphi_local_idx(0) >= m_num_nodes) {
                node_bucket_handle_y(0)++;
                sphi_local_idx(0) -= m_num_nodes;
              }

              if (sphi_local_idx(2) >= m_num_nodes) {
                node_bucket_handle_y(2)++;
                sphi_local_idx(2) -= m_num_nodes;
              }

              if (node_bucket_handle_y(0) < 0 ||
                  node_bucket_handle_y(0) >= m_particle_buckets.dim_size(0) ||
                  node_bucket_handle_y(1) < 0 ||
                  node_bucket_handle_y(1) >= m_particle_buckets.dim_size(1) ||
                  node_bucket_handle_y(2) < 0 ||
                  node_bucket_handle_y(2) >= m_particle_buckets.dim_size(2)) {
                bucket_node_idx_solid_phi_y(node_idx * 8 + (r * 2 + s) * 2 +
                                            0) = -1;
                bucket_node_idx_solid_phi_y(node_idx * 8 + (r * 2 + s) * 2 +
                                            1) = -1;
              } else {
                const int sphi_idx =
                    sphi_local_idx(2) * m_num_nodes * m_num_nodes +
                    sphi_local_idx(1) * m_num_nodes + sphi_local_idx(0);

                int nb_bucket_idx =
                    m_particle_buckets.bucket_index(node_bucket_handle_y);
                if (!m_bucket_activated[nb_bucket_idx]) {
                  bucket_node_idx_solid_phi_y(node_idx * 8 + (r * 2 + s) * 2 +
                                              0) = -1;
                  bucket_node_idx_solid_phi_y(node_idx * 8 + (r * 2 + s) * 2 +
                                              1) = -1;
                } else {
                  bucket_node_idx_solid_phi_y(node_idx * 8 + (r * 2 + s) * 2 +
                                              0) = nb_bucket_idx;
                  bucket_node_idx_solid_phi_y(node_idx * 8 + (r * 2 + s) * 2 +
                                              1) = sphi_idx;
                }
              }
            }

          for (int r = 0; r < 2; ++r)
            for (int s = 0; s < 2; ++s) {
              Vec3i node_bucket_handle_z = bucket_handle;

              Vec3i sphi_local_idx = Vec3i(i + r, j + s, k);

              if (sphi_local_idx(0) >= m_num_nodes) {
                node_bucket_handle_z(0)++;
                sphi_local_idx(0) -= m_num_nodes;
              }

              if (sphi_local_idx(1) >= m_num_nodes) {
                node_bucket_handle_z(1)++;
                sphi_local_idx(1) -= m_num_nodes;
              }

              if (node_bucket_handle_z(0) < 0 ||
                  node_bucket_handle_z(0) >= m_particle_buckets.dim_size(0) ||
                  node_bucket_handle_z(1) < 0 ||
                  node_bucket_handle_z(1) >= m_particle_buckets.dim_size(1) ||
                  node_bucket_handle_z(2) < 0 ||
                  node_bucket_handle_z(2) >= m_particle_buckets.dim_size(2)) {
                bucket_node_idx_solid_phi_z(node_idx * 8 + (r * 2 + s) * 2 +
                                            0) = -1;
                bucket_node_idx_solid_phi_z(node_idx * 8 + (r * 2 + s) * 2 +
                                            1) = -1;
              } else {
                const int sphi_idx =
                    sphi_local_idx(2) * m_num_nodes * m_num_nodes +
                    sphi_local_idx(1) * m_num_nodes + sphi_local_idx(0);

                int nb_bucket_idx =
                    m_particle_buckets.bucket_index(node_bucket_handle_z);
                if (!m_bucket_activated[nb_bucket_idx]) {
                  bucket_node_idx_solid_phi_z(node_idx * 8 + (r * 2 + s) * 2 +
                                              0) = -1;
                  bucket_node_idx_solid_phi_z(node_idx * 8 + (r * 2 + s) * 2 +
                                              1) = -1;
                } else {
                  bucket_node_idx_solid_phi_z(node_idx * 8 + (r * 2 + s) * 2 +
                                              0) = nb_bucket_idx;
                  bucket_node_idx_solid_phi_z(node_idx * 8 + (r * 2 + s) * 2 +
                                              1) = sphi_idx;
                }
              }
            }
        }
  });
}

void LiquidSimulator::connectPressureNodes() {
  const Vec3i ppdir[] = {Vec3i(-1, 0, 0),  Vec3i(1, 0, 0),   Vec3i(0, -1, 0),
                         Vec3i(0, 1, 0),   Vec3i(0, 0, -1),  Vec3i(0, 0, 1),
                         Vec3i(-1, -1, 0), Vec3i(1, -1, 0),  Vec3i(-1, 1, 0),
                         Vec3i(1, 1, 0),   Vec3i(-1, 0, -1), Vec3i(1, 0, -1),
                         Vec3i(0, -1, -1), Vec3i(0, 1, -1),  Vec3i(-1, 0, 1),
                         Vec3i(1, 0, 1),   Vec3i(0, -1, 1),  Vec3i(0, 1, 1)};

  m_particle_buckets.for_each_bucket_colored([&](int bucket_idx) {
    Vec3i bucket_handle = m_particle_buckets.bucket_handle(bucket_idx);
    if (!m_bucket_activated[bucket_idx]) return;

    VecXi& bucket_node_pressure_neighbors =
        m_node_pressure_neighbors[bucket_idx];
    VecXi& bucket_node_pp_neighbors = m_node_pp_neighbors[bucket_idx];

    const int count_p = getNumNodes(bucket_idx);
    bucket_node_pressure_neighbors.resize(count_p * 6 * 2);
    bucket_node_pp_neighbors.resize(count_p * 18 * 2);

    bucket_node_pressure_neighbors.setConstant(-1);
    bucket_node_pp_neighbors.setConstant(-1);

    for (int k = 0; k < m_num_nodes; ++k)
      for (int j = 0; j < m_num_nodes; ++j)
        for (int i = 0; i < m_num_nodes; ++i) {
          int node_idx = k * m_num_nodes * m_num_nodes + j * m_num_nodes + i;

          for (int r = 0; r < 18; ++r) {
            Vec3i node_bucket_handle_p = bucket_handle;
            Vec3i mac_local_idx = Vec3i(i, j, k) + ppdir[r];

            for (int s = 0; s < 3; ++s) {
              if (mac_local_idx(s) >= m_num_nodes) {
                node_bucket_handle_p(s)++;
                mac_local_idx(s) -= m_num_nodes;
              } else if (mac_local_idx(s) < 0) {
                node_bucket_handle_p(s)--;
                mac_local_idx(s) += m_num_nodes;
              }
            }

            if (!m_particle_buckets.has_bucket(node_bucket_handle_p)) continue;

            const int nb_bucket_idx =
                m_particle_buckets.bucket_index(node_bucket_handle_p);

            if (!m_bucket_activated[nb_bucket_idx]) continue;

            const int mac_idx = mac_local_idx(2) * m_num_nodes * m_num_nodes +
                                mac_local_idx(1) * m_num_nodes +
                                mac_local_idx(0);

            bucket_node_pp_neighbors[node_idx * 36 + r * 2 + 0] = nb_bucket_idx;
            bucket_node_pp_neighbors[node_idx * 36 + r * 2 + 1] = mac_idx;
          }

          for (int r = 0; r < 2; ++r) {
            Vec3i node_bucket_handle_x = bucket_handle;

            Vec3i mac_local_idx = Vec3i(i + r, j, k);

            if (mac_local_idx(0) >= m_num_nodes) {
              node_bucket_handle_x(0)++;
              mac_local_idx(0) -= m_num_nodes;
            }

            if (!m_particle_buckets.has_bucket(node_bucket_handle_x)) continue;

            int nb_bucket_idx =
                m_particle_buckets.bucket_index(node_bucket_handle_x);

            if (!m_bucket_activated[nb_bucket_idx]) continue;

            const int mac_idx = mac_local_idx(2) * m_num_nodes * m_num_nodes +
                                mac_local_idx(1) * m_num_nodes +
                                mac_local_idx(0);

            bucket_node_pressure_neighbors[node_idx * 12 + r * 2 + 0] =
                nb_bucket_idx;
            bucket_node_pressure_neighbors[node_idx * 12 + r * 2 + 1] = mac_idx;

            VecXi& nb_node_idxp = m_node_index_pressure_x[nb_bucket_idx];
            nb_node_idxp[mac_idx * 4 + (1 - r) * 2 + 0] = bucket_idx;
            nb_node_idxp[mac_idx * 4 + (1 - r) * 2 + 1] = node_idx;
          }

          for (int r = 0; r < 2; ++r) {
            Vec3i node_bucket_handle_y = bucket_handle;

            Vec3i mac_local_idx = Vec3i(i, j + r, k);

            if (mac_local_idx(1) >= m_num_nodes) {
              node_bucket_handle_y(1)++;
              mac_local_idx(1) -= m_num_nodes;
            }

            if (!m_particle_buckets.has_bucket(node_bucket_handle_y)) continue;

            int nb_bucket_idx =
                m_particle_buckets.bucket_index(node_bucket_handle_y);

            if (!m_bucket_activated[nb_bucket_idx]) continue;

            int mac_idx = mac_local_idx(2) * m_num_nodes * m_num_nodes +
                          mac_local_idx(1) * m_num_nodes + mac_local_idx(0);

            bucket_node_pressure_neighbors[node_idx * 12 + 4 + r * 2 + 0] =
                nb_bucket_idx;
            bucket_node_pressure_neighbors[node_idx * 12 + 4 + r * 2 + 1] =
                mac_idx;

            VecXi& nb_node_idxp = m_node_index_pressure_y[nb_bucket_idx];
            nb_node_idxp[mac_idx * 4 + (1 - r) * 2 + 0] = bucket_idx;
            nb_node_idxp[mac_idx * 4 + (1 - r) * 2 + 1] = node_idx;
          }

          for (int r = 0; r < 2; ++r) {
            Vec3i node_bucket_handle_z = bucket_handle;

            Vec3i mac_local_idx = Vec3i(i, j, k + r);

            if (mac_local_idx(2) >= m_num_nodes) {
              node_bucket_handle_z(2)++;
              mac_local_idx(2) -= m_num_nodes;
            }

            if (!m_particle_buckets.has_bucket(node_bucket_handle_z)) continue;

            int nb_bucket_idx =
                m_particle_buckets.bucket_index(node_bucket_handle_z);

            if (!m_bucket_activated[nb_bucket_idx]) continue;

            int mac_idx = mac_local_idx(2) * m_num_nodes * m_num_nodes +
                          mac_local_idx(1) * m_num_nodes + mac_local_idx(0);

            bucket_node_pressure_neighbors[node_idx * 12 + 8 + r * 2 + 0] =
                nb_bucket_idx;
            bucket_node_pressure_neighbors[node_idx * 12 + 8 + r * 2 + 1] =
                mac_idx;

            VecXi& nb_node_idxp = m_node_index_pressure_z[nb_bucket_idx];
            nb_node_idxp[mac_idx * 4 + (1 - r) * 2 + 0] = bucket_idx;
            nb_node_idxp[mac_idx * 4 + (1 - r) * 2 + 1] = node_idx;
          }
        }
  });
}

void LiquidSimulator::postAllocateNodes() {
  const int num_buckets = m_particle_buckets.size();

  allocateNodeVectors(m_node_mass_fluid_x, m_node_mass_fluid_y,
                      m_node_mass_fluid_z);
  allocateNodeVectors(m_node_vel_fluid_x, m_node_vel_fluid_y,
                      m_node_vel_fluid_z);
  allocateNodeVectors(m_node_vel_fluid_plus_x, m_node_vel_fluid_plus_y,
                      m_node_vel_fluid_plus_z);
  allocateNodeVectors(m_node_vol_fluid_x, m_node_vol_fluid_y,
                      m_node_vol_fluid_z);
  allocateNodeVectors(m_node_rhs_fluid_x, m_node_rhs_fluid_y,
                      m_node_rhs_fluid_z);
  allocateNodeVectors(m_node_vel_elastic_x, m_node_vel_elastic_y,
                      m_node_vel_elastic_z);
  allocateNodeVectors(m_node_drag_coeff_x, m_node_drag_coeff_y,
                      m_node_drag_coeff_z);
  allocateNodeVectors(m_node_vel_constraint_coeff_x,
                      m_node_vel_constraint_coeff_y,
                      m_node_vel_constraint_coeff_z);
  allocateNodeVectors(m_node_pressure_grad_x, m_node_pressure_grad_y,
                      m_node_pressure_grad_z);
  allocateNodeVectors(m_node_solid_vel_x, m_node_solid_vel_y,
                      m_node_solid_vel_z);
  allocateNodeVectors(m_node_liquid_valid_x, m_node_liquid_valid_y,
                      m_node_liquid_valid_z);

  allocateNodeVectors(m_node_vol_change_p);
  allocateNodeVectors(m_node_elastic_vf_p);
  allocateNodeVectors(m_node_extra_flow_vol_p);
  allocateNodeVectors(m_node_solid_phi);
  allocateNodeVectors(m_node_pressure);

  allocateNodeVectors(m_node_components, m_liquid_info.num_components);

  // only needs 6 due to symmetry
  //        allocateNodeVectors(m_node_pressure_hessian_p, 6);
}

void LiquidSimulator::resampleNodes() {
  const int num_parts = numParticles();

  preAllocateNodes();

  if (num_parts > 0) {
    auto part_pos = [this](int pidx) -> Vec3x {
      return m_x.segment<3>(pidx * 3);
    };

    findNodes(m_particle_buckets, m_particle_nodes_x, Vec3x(0.0, 0.5, 0.5),
              part_pos);
    findNodes(m_particle_buckets, m_particle_nodes_y, Vec3x(0.5, 0.0, 0.5),
              part_pos);
    findNodes(m_particle_buckets, m_particle_nodes_z, Vec3x(0.5, 0.5, 0.0),
              part_pos);
    findNodes(m_particle_buckets, m_particle_nodes_p, Vec3x(0.5, 0.5, 0.5),
              part_pos);
  }

  if (m_strands_global_local.size() > 0) {
    auto vert_pos = [this](int pidx) -> Vec3x {
      const pair<int, int>& local_idx = m_strands_global_local[pidx];
      const ElasticStrand* strand = m_strands[local_idx.first];
      return strand->getVertex(local_idx.second);
    };

    findNodes(m_vertex_buckets, m_vertex_nodes_x, Vec3x(0.0, 0.5, 0.5),
              vert_pos);
    findNodes(m_vertex_buckets, m_vertex_nodes_y, Vec3x(0.5, 0.0, 0.5),
              vert_pos);
    findNodes(m_vertex_buckets, m_vertex_nodes_z, Vec3x(0.5, 0.5, 0.0),
              vert_pos);
    findNodes(m_vertex_buckets, m_vertex_nodes_p, Vec3x(0.5, 0.5, 0.5),
              vert_pos);
  }

  expandFluidNodesMarked(1);

  // generate nodes in all activated buckets
  generateNodes();

  connectSolidPhiNodes();
  // connect pressure node to MAC nodes
  connectPressureNodes();

  postAllocateNodes();
}

bool LiquidSimulator::computeStress() {
  // update Fe, J, b
  return true;
}

void LiquidSimulator::updateVertexWeights() {
  const Scalar h = getDx();

  const int num_vertex = m_strands_global_local.size();

  if (!num_vertex) return;

  auto accum_weight = [&](const Mat27x2i& indices, const Vec3x& pos,
                          Mat27x4f& weights, int ir) -> bool {
    bool touched = false;

    for (int nidx = 0; nidx < indices.rows(); ++nidx) {
      const int node_bucket_idx = indices(nidx, 0);
      const int node_idx = indices(nidx, 1);

      Scalar w = 0.0;

      if (node_bucket_idx != -1 && node_idx != -1 &&
          m_bucket_activated[node_bucket_idx]) {
        const Vec3x& np = getNodePos(node_bucket_idx, node_idx, ir);
        Vec3x dx = (pos - np) / h;
        w = N_kernel<2>(dx);
      }

      weights(nidx, ir) = w;

      if (w > 0.0) {
        touched = true;
      }
    }

    return touched;
  };

  vector<unsigned char> vertex_touched(num_vertex);

  for_each(0, num_vertex, [&](int pidx) {
    const Mat27x2i& indices_x = m_vertex_nodes_x[pidx];
    const Mat27x2i& indices_y = m_vertex_nodes_y[pidx];
    const Mat27x2i& indices_z = m_vertex_nodes_z[pidx];
    const Mat27x2i& indices_p = m_vertex_nodes_p[pidx];

    Mat27x4f& weights = m_vertex_weights[pidx];

    const pair<int, int>& local_idx = m_strands_global_local[pidx];
    const ElasticStrand* strand = m_strands[local_idx.first];
    const Vec3x pos = strand->getVertex(local_idx.second);

    unsigned char touched = false;

    touched |= accum_weight(indices_x, pos, weights, 0);
    touched |= accum_weight(indices_y, pos, weights, 1);
    touched |= accum_weight(indices_z, pos, weights, 2);
    touched |= accum_weight(indices_p, pos, weights, 3);

    vertex_touched[pidx] = touched;
  });

  //        for(int i = 0; i < num_vertex; ++i) {
  //            std::cout << "vt " << i << ": " << (int) vertex_touched[i] <<
  //            std::endl;
  //        }

  const int num_strands = (int)m_strands.size();
  for_each(0, num_strands, [&](int strand_idx) {
    const int global_base_idx = m_strands_local_global_base[strand_idx];
    const int num_verts = m_strands[strand_idx]->getNumVertices();

    unsigned char touched = false;
    for (int i = 0; i < num_verts && !touched; ++i) {
      touched |= vertex_touched[global_base_idx + i];
    }
    m_strands_submerged[strand_idx] = touched;
  });
  //
  //        for(int i = 0; i < num_strands; ++i) {
  //            std::cout << "sbm " << i << ": " << (int) m_strands_submerged[i]
  //            << std::endl;
  //        }
}

void LiquidSimulator::findInterfacingSegments() {
  const int num_strand = m_strands.size();
  if (!num_strand) return;

  std::vector<std::vector<std::pair<int, int> > > strand_inside_edges(
      num_strand);

  for_each(0, num_strand, [&](int strand_idx) {
    ElasticStrand* strand = m_strands[strand_idx];
    strand->dynamics().updateInterfaceSegments();

    const std::vector<int>& inside_edges = strand->dynamics().getInsideEdges();
    const int num_edges = inside_edges.size();

    strand_inside_edges[strand_idx].resize(num_edges);

    for (int i = 0; i < num_edges; ++i) {
      strand_inside_edges[strand_idx][i] =
          std::pair<int, int>(strand_idx, inside_edges[i]);
    }
  });

  m_strands_inside_segs.resize(0);
  parallel_concatenate(m_strands_inside_segs, strand_inside_edges);

  const int num_inside_edges = m_strands_inside_segs.size();

  auto edge_pos = [this](int eidx) -> Vec3x {
    const pair<int, int>& local_idx = m_strands_inside_segs[eidx];
    const ElasticStrand* strand = m_strands[local_idx.first];
    return (strand->getVertex(local_idx.second) +
            strand->getVertex(local_idx.second + 1)) *
           0.5;
  };

  m_edge_buckets.resize(m_particle_buckets.ni, m_particle_buckets.nj,
                        m_particle_buckets.nk);
  m_edge_buckets.sort(num_inside_edges, [&](int eidx, int& i, int& j, int& k) {
    const Vec3x pos = edge_pos(eidx);
    i = (int)std::floor((pos(0) - m_grid_mincorner(0)) / m_bucket_size);
    j = (int)std::floor((pos(1) - m_grid_mincorner(1)) / m_bucket_size);
    k = (int)std::floor((pos(2) - m_grid_mincorner(2)) / m_bucket_size);
  });

  m_inside_edge_nodes_x.resize(num_inside_edges);
  m_inside_edge_nodes_y.resize(num_inside_edges);
  m_inside_edge_nodes_z.resize(num_inside_edges);

  findNodes(m_edge_buckets, m_inside_edge_nodes_x, Vec3x(0.0, 0.5, 0.5),
            edge_pos, false);
  findNodes(m_edge_buckets, m_inside_edge_nodes_y, Vec3x(0.5, 0.0, 0.5),
            edge_pos, false);
  findNodes(m_edge_buckets, m_inside_edge_nodes_z, Vec3x(0.5, 0.5, 0.0),
            edge_pos, false);

  updateEdgeWeights();
  buildNodeEdgePairs();
}

void LiquidSimulator::updateEdgeWeights() {
  const Scalar h = getDx();

  auto accum_weight = [&](const Mat27x2i& indices, const Vec3x& pos,
                          Mat27x3f& weights, int ir) {
    for (int nidx = 0; nidx < indices.rows(); ++nidx) {
      const int node_bucket_idx = indices(nidx, 0);
      const int node_idx = indices(nidx, 1);

      Scalar w = 0.0;

      if (node_bucket_idx != -1 && node_idx != -1 &&
          m_bucket_activated[node_bucket_idx]) {
        const Vec3x& np = getNodePos(node_bucket_idx, node_idx, ir);
        Vec3x dx = (pos - np) / h;
        w = N_kernel<2>(dx);
      }

      weights(nidx, ir) = w;
    }
  };

  auto edge_pos = [this](int eidx) -> Vec3x {
    const pair<int, int>& local_idx = m_strands_inside_segs[eidx];
    const ElasticStrand* strand = m_strands[local_idx.first];
    return (strand->getVertex(local_idx.second) +
            strand->getVertex(local_idx.second + 1)) *
           0.5;
  };

  const int num_edges = m_strands_inside_segs.size();

  m_inside_edge_weights.resize(num_edges);

  for_each(0, num_edges, [&](int pidx) {
    const Mat27x2i& indices_x = m_inside_edge_nodes_x[pidx];
    const Mat27x2i& indices_y = m_inside_edge_nodes_y[pidx];
    const Mat27x2i& indices_z = m_inside_edge_nodes_z[pidx];

    Mat27x3f& weights = m_inside_edge_weights[pidx];

    const Vec3x pos = edge_pos(pidx);

    accum_weight(indices_x, pos, weights, 0);
    accum_weight(indices_y, pos, weights, 1);
    accum_weight(indices_z, pos, weights, 2);
  });
}

void LiquidSimulator::updateParticleWeights() {
  const Scalar h = getDx();

  auto accum_weight = [&](const Mat27x2i& indices, const Vec3x& pos,
                          Mat27x4f& weights, int ir) {
    for (int nidx = 0; nidx < indices.rows(); ++nidx) {
      const int node_bucket_idx = indices(nidx, 0);
      const int node_idx = indices(nidx, 1);

      const Vec3x& np = getNodePos(node_bucket_idx, node_idx, ir);
      Vec3x dx = (pos - np) / h;

      weights(nidx, ir) = N_kernel<2>(dx);
    }
  };

  for_each(0, numParticles(), [&](int pidx) {
    const Mat27x2i& indices_x = m_particle_nodes_x[pidx];
    const Mat27x2i& indices_y = m_particle_nodes_y[pidx];
    const Mat27x2i& indices_z = m_particle_nodes_z[pidx];
    const Mat27x2i& indices_p = m_particle_nodes_p[pidx];

    Mat27x4f& weights = m_particle_weights[pidx];
    //            Mat27x3f& grad_x = m_particle_grad_x[pidx];
    //            Mat27x3f& grad_y = m_particle_grad_y[pidx];
    //            Mat27x3f& grad_z = m_particle_grad_z[pidx];

    const Vec3x& pos = m_x.segment<3>(pidx * 3);

    //            accum_weight_grad(indices_x, pos, weights, grad_x, 0);
    //            accum_weight_grad(indices_y, pos, weights, grad_y, 1);
    //            accum_weight_grad(indices_z, pos, weights, grad_z, 2);
    accum_weight(indices_x, pos, weights, 0);
    accum_weight(indices_y, pos, weights, 1);
    accum_weight(indices_z, pos, weights, 2);
    accum_weight(indices_p, pos, weights, 3);
  });
}

template <int N>
void LiquidSimulator::buildNodePairs(
    const Sorter& buckets, NodeObjIndex& node_obj,
    const vector<Eigen::Matrix<float, 27, N> >& weights,
    const vector<Mat27x2i>& obj_node, int isel, std::function<bool(int)> func) {
  const int num_buckets = (int)buckets.size();

  if ((int)node_obj.size() != num_buckets) node_obj.resize(num_buckets);

  // re-allocate space
  buckets.for_each_bucket([&](int bucket_idx) {
    int num_nodes = getNumNodes(bucket_idx);

    auto& bucket_node_particles = node_obj[bucket_idx];
    if ((int)bucket_node_particles.size() != num_nodes)
      bucket_node_particles.resize(num_nodes);

    for (int i = 0; i < num_nodes; ++i) {
      bucket_node_particles[i].resize(0);
    }
  });

  if (!weights.size()) return;

  buckets.for_each_bucket_particles_colored(
      [&](int pidx, int bucket_idx) {
        if (!m_bucket_activated[bucket_idx]) return;
        if (func && !func(pidx)) return;

        auto& indices = obj_node[pidx];

        auto& w = weights[pidx];

        for (int i = 0; i < indices.rows(); ++i) {
          if (indices(i, 0) >= 0 && indices(i, 1) >= 0 &&
              m_bucket_activated[indices(i, 0)] && w(i, isel) > 0.0) {
            node_obj[indices(i, 0)][indices(i, 1)].emplace_back(
                std::pair<int, short>(pidx, i));
          }
        }
      },
      3);
}

void LiquidSimulator::buildNodeEdgePairs() {
  auto func = [this](int pidx) -> bool {
    const int strand_idx = m_strands_inside_segs[pidx].first;
    const int edge_idx = m_strands_inside_segs[pidx].second;

    return m_strands[strand_idx]->dynamics().isEdgeInterfacing(edge_idx);
  };

  buildNodePairs<3>(m_edge_buckets, m_node_edge_x, m_inside_edge_weights,
                    m_inside_edge_nodes_x, 0, func);
  buildNodePairs<3>(m_edge_buckets, m_node_edge_y, m_inside_edge_weights,
                    m_inside_edge_nodes_y, 1, func);
  buildNodePairs<3>(m_edge_buckets, m_node_edge_z, m_inside_edge_weights,
                    m_inside_edge_nodes_z, 2, func);
}

void LiquidSimulator::buildNodeVertexPairs() {
  buildNodePairs<4>(m_vertex_buckets, m_node_vertex_x, m_vertex_weights,
                    m_vertex_nodes_x, 0);
  buildNodePairs<4>(m_vertex_buckets, m_node_vertex_y, m_vertex_weights,
                    m_vertex_nodes_y, 1);
  buildNodePairs<4>(m_vertex_buckets, m_node_vertex_z, m_vertex_weights,
                    m_vertex_nodes_z, 2);
  buildNodePairs<4>(m_vertex_buckets, m_node_vertex_p, m_vertex_weights,
                    m_vertex_nodes_p, 3);
}

void LiquidSimulator::buildNodeParticlePairs() {
  // std::cout << "[buildNodeParticlePairs 0]" << std::endl;
  buildNodePairs<4>(m_particle_buckets, m_node_particles_x, m_particle_weights,
                    m_particle_nodes_x, 0);

  // std::cout << "[buildNodeParticlePairs 1]" << std::endl;
  buildNodePairs<4>(m_particle_buckets, m_node_particles_y, m_particle_weights,
                    m_particle_nodes_y, 1);

  // std::cout << "[buildNodeParticlePairs 2]" << std::endl;
  buildNodePairs<4>(m_particle_buckets, m_node_particles_z, m_particle_weights,
                    m_particle_nodes_z, 2);

  // std::cout << "[buildNodeParticlePairs 3]" << std::endl;
  buildNodePairs<4>(m_particle_buckets, m_node_particles_p, m_particle_weights,
                    m_particle_nodes_p, 3);

  // std::cout << "[buildNodeParticlePairs 4]" << std::endl;
}

bool LiquidSimulator::updatePlasticity() {
  if (!m_liquid_info.solve_viscosity) return true;

  for_each(0, numParticles(), [&](int pidx) {
    const Mat3x& F = m_Fe.block<3, 3>(pidx * 3, 0);
    // check_isnan("F in UP", F.sum());

    const Mat3x& be_bar = m_b.block<3, 3>(pidx * 3, 0);
    // check_isnan("bebar in UP", be_bar.sum());

    const Scalar& J = m_J(pidx);
    // check_isnan("J in UP", J);

    const Mat3x inv_be_bar = be_bar.inverse();
    // check_isnan("inv_bebar in UP", inv_be_bar.sum());

    const Mat3x Cp = F.transpose() * inv_be_bar * F / pow(J, 2.0 / 3.0);

    // check_isnan("Cp in UP", Cp.sum());

    Eigen::SelfAdjointEigenSolver<Mat3x> es;
    es.compute(Cp);  // Cp = U P U^T
    const Mat3x U = es.eigenvectors();
    Mat3x P = es.eigenvalues().asDiagonal();

    P(0, 0) = std::max(0.01, P(0, 0));
    P(1, 1) = std::max(0.01, P(1, 1));
    P(2, 2) = std::max(0.01, P(2, 2));

    // check_isnan("P in UP", P.sum());
    // check_isnan("U in UP", U.sum());

    Mat3x newP = Mat3x::Zero();

    const VecXx& color = m_components.segment(
        pidx * m_liquid_info.num_components, m_liquid_info.num_components);
    // check_isnan("color in UP", color);

    const Scalar plastic_relax_coeff =
        m_liquid_info.plastic_relax_coeff.dot(color);
    const Scalar P_coeff = (exp(-m_dt / plastic_relax_coeff) - 1.0) * 0.5;

    newP(0, 0) = pow(P(0, 0), P_coeff);
    newP(1, 1) = pow(P(1, 1), P_coeff);
    newP(2, 2) = pow(P(2, 2), P_coeff);

    const Mat3x newF = F * U * newP * U.transpose();
    //            if(std::isnan(newF.sum())){
    //                std::cout << color << std::endl << std::endl;
    //                std::cout << plastic_relax_coeff << std::endl;
    //                std::cout << P_coeff << std::endl << std::endl;
    //                std::cout << P << std::endl << std::endl;
    //                std::cout << U << std::endl << std::endl;
    //                std::cout << F << std::endl << std::endl;
    //            }
    //
    // check_isnan("newF in UP", newF.sum());

    const Scalar newJ = newF.determinant();
    if (newJ == 0.0 || std::isnan(newJ)) {
#pragma omp critical
      {
        std::cout << "color: " << color << std::endl << std::endl;
        std::cout << "relaxcoeff: " << plastic_relax_coeff << std::endl
                  << std::endl;
        std::cout << "Pcoeff: " << P_coeff << std::endl << std::endl;
        std::cout << "P: " << P << std::endl << std::endl;
        std::cout << "newP: " << newP << std::endl << std::endl;
        std::cout << "U: " << U << std::endl << std::endl;
        std::cout << "F: " << F << std::endl << std::endl;
        std::cout << "newF: " << newF << std::endl << std::endl;

        std::cout << "J is bad in UP! " << newJ << std::endl;
        exit(-1);
      }
    }

    // check_isnan("newJ in UP", newJ);

    m_Fe.block<3, 3>(pidx * 3, 0) = pow(J / newJ, 1.0 / 3.0) * newF;
    // check_isnan("m_Fe in UP", m_Fe.block<3, 3>(pidx * 3, 0).sum());
  });
  return true;
}

bool LiquidSimulator::updateSolidPhi() {
  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    const int num_nodes = getNumNodes(bucket_idx);

    VecXx& node_phi = m_node_solid_phi[bucket_idx];

    VecXx& node_solid_vel_x = m_node_solid_vel_x[bucket_idx];
    VecXx& node_solid_vel_y = m_node_solid_vel_y[bucket_idx];
    VecXx& node_solid_vel_z = m_node_solid_vel_z[bucket_idx];

    for (int i = 0; i < num_nodes; ++i) {
      Vec3x vel;
      node_phi(i) = DistanceFieldObject::computePhi(
          getNodePosSolidPhi(bucket_idx, i), m_fields);

      DistanceFieldObject::computePhiVel(getNodePosX(bucket_idx, i), vel,
                                         m_fields);
      node_solid_vel_x(i) = vel(0);

      DistanceFieldObject::computePhiVel(getNodePosY(bucket_idx, i), vel,
                                         m_fields);
      node_solid_vel_y(i) = vel(1);

      DistanceFieldObject::computePhiVel(getNodePosZ(bucket_idx, i), vel,
                                         m_fields);
      node_solid_vel_z(i) = vel(2);
    }
  });

  return true;
}

void LiquidSimulator::updateSolidWeights() {
  const int num_buckets = m_particle_buckets.size();
  m_node_liquid_weight_x.resize(num_buckets);
  m_node_liquid_weight_y.resize(num_buckets);
  m_node_liquid_weight_z.resize(num_buckets);

  const Scalar dx = getDx();

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    const VecXi& bucket_node_idx_solid_phi_x =
        m_node_index_solid_phi_x[bucket_idx];
    const VecXi& bucket_node_idx_solid_phi_y =
        m_node_index_solid_phi_y[bucket_idx];
    const VecXi& bucket_node_idx_solid_phi_z =
        m_node_index_solid_phi_z[bucket_idx];

    const int num_solid_phi = getNumNodes(bucket_idx);

    VecXx& bucket_weight_x = m_node_liquid_weight_x[bucket_idx];
    VecXx& bucket_weight_y = m_node_liquid_weight_y[bucket_idx];
    VecXx& bucket_weight_z = m_node_liquid_weight_z[bucket_idx];

    if (bucket_weight_x.size() != num_solid_phi)
      bucket_weight_x.resize(num_solid_phi);
    if (bucket_weight_y.size() != num_solid_phi)
      bucket_weight_y.resize(num_solid_phi);
    if (bucket_weight_z.size() != num_solid_phi)
      bucket_weight_z.resize(num_solid_phi);

    for (int i = 0; i < num_solid_phi; ++i) {
      const Vec8i& indices = bucket_node_idx_solid_phi_x.segment<8>(i * 8);
      Scalar phi0 = 0.5 * dx;
      Scalar phi1 = 0.5 * dx;
      Scalar phi2 = 0.5 * dx;
      Scalar phi3 = 0.5 * dx;

      if (indices[0] >= 0 && m_bucket_activated[indices[0]])
        phi0 = m_node_solid_phi[indices[0]][indices[1]];
      if (indices[2] >= 0 && m_bucket_activated[indices[2]])
        phi1 = m_node_solid_phi[indices[2]][indices[3]];
      if (indices[4] >= 0 && m_bucket_activated[indices[4]])
        phi2 = m_node_solid_phi[indices[4]][indices[5]];
      if (indices[6] >= 0 && m_bucket_activated[indices[6]])
        phi3 = m_node_solid_phi[indices[6]][indices[7]];

      bucket_weight_x(i) =
          clamp(1.0 - fraction_inside(phi0, phi1, phi2, phi3), 0.0, 1.0);
    }

    for (int i = 0; i < num_solid_phi; ++i) {
      const Vec8i& indices = bucket_node_idx_solid_phi_y.segment<8>(i * 8);
      Scalar phi0 = 0.5 * dx;
      Scalar phi1 = 0.5 * dx;
      Scalar phi2 = 0.5 * dx;
      Scalar phi3 = 0.5 * dx;

      if (indices[0] >= 0 && m_bucket_activated[indices[0]])
        phi0 = m_node_solid_phi[indices[0]][indices[1]];
      if (indices[2] >= 0 && m_bucket_activated[indices[2]])
        phi1 = m_node_solid_phi[indices[2]][indices[3]];
      if (indices[4] >= 0 && m_bucket_activated[indices[4]])
        phi2 = m_node_solid_phi[indices[4]][indices[5]];
      if (indices[6] >= 0 && m_bucket_activated[indices[6]])
        phi3 = m_node_solid_phi[indices[6]][indices[7]];

      bucket_weight_y(i) =
          clamp(1.0 - fraction_inside(phi0, phi1, phi2, phi3), 0.0, 1.0);
    }

    for (int i = 0; i < num_solid_phi; ++i) {
      const Vec8i& indices = bucket_node_idx_solid_phi_z.segment<8>(i * 8);
      Scalar phi0 = 0.5 * dx;
      Scalar phi1 = 0.5 * dx;
      Scalar phi2 = 0.5 * dx;
      Scalar phi3 = 0.5 * dx;

      if (indices[0] >= 0 && m_bucket_activated[indices[0]])
        phi0 = m_node_solid_phi[indices[0]][indices[1]];
      if (indices[2] >= 0 && m_bucket_activated[indices[2]])
        phi1 = m_node_solid_phi[indices[2]][indices[3]];
      if (indices[4] >= 0 && m_bucket_activated[indices[4]])
        phi2 = m_node_solid_phi[indices[4]][indices[5]];
      if (indices[6] >= 0 && m_bucket_activated[indices[6]])
        phi3 = m_node_solid_phi[indices[6]][indices[7]];

      bucket_weight_z(i) =
          clamp(1.0 - fraction_inside(phi0, phi1, phi2, phi3), 0.0, 1.0);
    }
  });
}

void LiquidSimulator::updateLiquidPhi() {
  const int num_buckets = (int)m_particle_buckets.size();

  m_node_liquid_phi.resize(num_buckets);

  const Scalar dx = getDx();

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    const int num_nodes = getNumNodes(bucket_idx);

    m_node_liquid_phi[bucket_idx].resize(num_nodes);
    m_node_liquid_phi[bucket_idx].setConstant(3.0 * m_bucket_size);
  });

  if (numParticles() == 0) return;

  m_particle_buckets.for_each_bucket_particles_colored(
      [&](int pidx, int) {
        const auto& indices = m_particle_nodes_p[pidx];

        const Vec3x& pos =
            m_x.segment<3>(pidx * 3);  // + m_fluid_v.segment<3>(pidx * 4) * dt;

        for (int i = 0; i < indices.rows(); ++i) {
          if (!m_bucket_activated[indices(i, 0)]) continue;

          VecXx& phis = m_node_liquid_phi[indices(i, 0)];
          assert(indices(i, 1) >= 0 && indices(i, 1) < phis.size());

          const Vec3x& np = getNodePosP(indices(i, 0), indices(i, 1));

          const Scalar phi =
              (pos - np).norm() - std::max(dx * 0.883644, m_radius(pidx));

          if (phi < phis(indices(i, 1))) {
            phis(indices(i, 1)) = phi;
          }
        }
      },
      3);

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    VecXx& bucket_liquid_phi = m_node_liquid_phi[bucket_idx];

    const int num_pressure = bucket_liquid_phi.size();

    for (int i = 0; i < num_pressure; ++i) {
      const Vec3x& np = getNodePosP(bucket_idx, i);

      const Scalar sphi = DistanceFieldObject::computePhi(np, m_fields);
      if (sphi < 0.0) bucket_liquid_phi(i) = -0.5 * dx;
    }
  });
}

void LiquidSimulator::updateOptiVolume() {
  const Scalar rad_fine = defaultRadiusMultiplier() * getDx() *
                          m_liquid_info.particle_cell_multiplier;
  const Scalar V_fine = 4.0 / 3.0 * M_PI * rad_fine * rad_fine * rad_fine;

  const int num_parts = numParticles();
  for_each(0, num_parts, [&](int pidx) {
    const Scalar mrel = m_vol(pidx) / V_fine;

    if (mrel < 0.5)
      m_classifier[pidx] = PC_S;
    else if (mrel <= 0.9)
      m_classifier[pidx] = PC_s;
    else if (mrel <= 1.1)
      m_classifier[pidx] = PC_o;
    else if (mrel <= 2.0)
      m_classifier[pidx] = PC_l;
    else
      m_classifier[pidx] = PC_L;
  });
}

void LiquidSimulator::particleCheckWeakened() {
  const int num_parts = numParticles();

  for_each(0, num_parts, [&](int pidx) {
    const Mat3x inv_be_bar = m_b.block<3, 3>(pidx * 3, 0).inverse();
    const Scalar& J = m_J(pidx);
    const Mat3x& F = m_Fe.block<3, 3>(pidx * 3, 0);
    const Mat3x Cp = F.transpose() * inv_be_bar * F * pow(J, -2.0 / 3.0);
    const Mat3x dev_Cp = Cp - Mat3x::Identity() * Cp.trace() / 3.0;
    const Scalar len_dev_Cp = dev_Cp.norm();
    const VecXx& color = m_components.segment(
        pidx * m_liquid_info.num_components, m_liquid_info.num_components);
    m_weakened[pidx] =
        len_dev_Cp >= m_liquid_info.plastic_weaken_strain.dot(color);
  });
}

void LiquidSimulator::splitLiquidParticles() {
  const int num_fluids = numParticles();
  if (!num_fluids) return;

  std::vector<std::vector<Vec3x> > new_part_pos(num_fluids);
  std::vector<int> n_additional(num_fluids, 0);

  const Scalar rad_fine = defaultRadiusMultiplier() * getDx() *
                          m_liquid_info.particle_cell_multiplier;
  const Scalar V_fine = 4.0 / 3.0 * M_PI * rad_fine * rad_fine * rad_fine;

  LockGuard lock(m_particle_mutex);

  for_each(0, num_fluids, [&](int pidx) {
    if (m_classifier[pidx] != PC_L) return;

    const int n_split = std::min((int)std::ceil(m_vol(pidx) / V_fine),
                                 (int)sphere_pattern::max_vector_length);
    if (n_split <= 1) return;

    const Vec3x center = m_x.segment<3>(pidx * 3);
    const Scalar rad = m_radius(pidx);

    const Scalar new_vol = m_vol(pidx) / (Scalar)n_split;
    const Scalar new_rest_vol = m_rest_vol(pidx) / (Scalar)n_split;
    const Scalar new_rad = pow(new_vol / M_PI * 0.75, 1.0 / 3.0);
    const Scalar splat_rad = std::max(new_rad, rad - new_rad) * 0.75;

    new_part_pos[pidx].resize(n_split - 1);

    Mat3x M = Mat3x::Random();
    Mat3x Q, R;
    QRDecompose<Scalar, 3>(M, Q, R);

    for (int i = 1; i < n_split; ++i) {
      new_part_pos[pidx][i - 1] =
          center + Q * m_sphere_pattern[n_split].segment<3>(i * 3) * splat_rad;
    }

    n_additional[pidx] = n_split - 1;

    m_x.segment<3>(pidx * 3) =
        center + Q * m_sphere_pattern[n_split].segment<3>(0) * splat_rad;
    m_radius(pidx) = new_rad;
    m_vol(pidx) = new_vol;
    m_rest_vol(pidx) = new_rest_vol;

    const VecXx& color = m_components.segment(
        pidx * m_liquid_info.num_components, m_liquid_info.num_components);
    m_m.segment<3>(pidx * 3).setConstant(new_vol * getDensity(color));
    m_classifier[pidx] = PC_o;
  });

  std::partial_sum(n_additional.begin(), n_additional.end(),
                   n_additional.begin());

  const int num_add = n_additional[n_additional.size() - 1];

  if (num_add == 0) return;

  const int old_num_parts = numParticles();

  conservativeResizeParticles(old_num_parts + num_add);

  for_each(0, num_fluids, [&](int pidx_parent) {
    const int idx_np = ((pidx_parent == 0) ? 0 : n_additional[pidx_parent - 1]);
    const int pidx_new_parts = idx_np + old_num_parts;
    const int num_new_parts = new_part_pos[pidx_parent].size();

    for (int i = 0; i < num_new_parts; ++i) {
      const int pidx = pidx_new_parts + i;
      m_x.segment<3>(pidx * 3) = new_part_pos[pidx_parent][i];
      m_v.segment<3>(pidx * 3) = m_v.segment<3>(pidx_parent * 3);
      m_m.segment<3>(pidx * 3) = m_m.segment<3>(pidx_parent * 3);
      m_vol(pidx) = m_vol(pidx_parent);
      m_rest_vol(pidx) = m_rest_vol(pidx_parent);
      m_radius(pidx) = m_radius(pidx_parent);
      m_particle_group[pidx] = m_particle_group[pidx_parent];
      m_B.block<3, 3>(pidx * 3, 0) = m_B.block<3, 3>(pidx_parent * 3, 0);
      m_b.block<3, 3>(pidx * 3, 0) = m_b.block<3, 3>(pidx_parent * 3, 0);
      m_b_plus.block<3, 3>(pidx * 3, 0) =
          m_b_plus.block<3, 3>(pidx_parent * 3, 0);
      m_Fe.block<3, 3>(pidx * 3, 0) = m_Fe.block<3, 3>(pidx_parent * 3, 0);
      m_Fe_plus.block<3, 3>(pidx * 3, 0) =
          m_Fe_plus.block<3, 3>(pidx_parent * 3, 0);
      m_b_trial.block<3, 3>(pidx * 3, 0) =
          m_b_trial.block<3, 3>(pidx_parent * 3, 0);
      m_proj_func.segment<3>(pidx * 3) =
          m_proj_func.segment<3>(pidx_parent * 3);
      m_J(pidx) = m_J(pidx_parent);
      m_weakened(pidx) = m_weakened(pidx_parent);
      m_classifier[pidx] = m_classifier[pidx_parent];
      m_components.segment(pidx * m_liquid_info.num_components,
                           m_liquid_info.num_components) =
          m_components.segment(pidx_parent * m_liquid_info.num_components,
                               m_liquid_info.num_components);
    }
  });

  updateOptiVolume();
}
void LiquidSimulator::mergeLiquidParticles() {
  if (!m_liquid_info.correction_step) return;

  const int num_parts = numParticles();
  if (!num_parts) return;

  std::vector<unsigned char> removed(num_parts, false);

  std::vector<Scalar> gathered_vol(num_parts, 0.0);
  std::vector<Vec3x> gathered_moment(num_parts, Vec3x::Zero());
  std::vector<Vec3x> gathered_grad_p(num_parts, Vec3x::Zero());
  std::vector<Mat3x> gathered_L(num_parts, Mat3x::Zero());
  std::vector<Mat3x> gathered_def_grad(num_parts, Mat3x::Zero());
  std::vector<Mat3x> gathered_b(num_parts, Mat3x::Zero());
  std::vector<VecXx> gathered_color(num_parts,
                                    VecXx::Zero(m_liquid_info.num_components));

  const Scalar rad_fine = defaultRadiusMultiplier() * getDx() *
                          m_liquid_info.particle_cell_multiplier;
  const Scalar V_fine = 4.0 / 3.0 * M_PI * rad_fine * rad_fine * rad_fine;

  const int correction_selector = rand() % m_liquid_info.correction_step;

  // prepare sorter
  m_particle_buckets.sort(
      numParticles(), [&](int pidx, int& i, int& j, int& k) {
        i = (int)std::floor((m_x(pidx * 3 + 0) - m_grid_mincorner(0)) /
                            m_bucket_size);
        j = (int)std::floor((m_x(pidx * 3 + 1) - m_grid_mincorner(1)) /
                            m_bucket_size);
        k = (int)std::floor((m_x(pidx * 3 + 2) - m_grid_mincorner(2)) /
                            m_bucket_size);
      });

  m_particle_buckets.for_each_bucket_particles_colored_randomized(
      [&](int pidx, int bucket_idx) {
        if (m_classifier[pidx] != PC_S && m_classifier[pidx] != PC_l) return;

        if (pidx % m_liquid_info.correction_step != correction_selector) return;

        Scalar should_rad = rad_fine * 2.0;

        if (m_classifier[pidx] == PC_S) {
          // try upgrade level
          const Scalar full_vol = m_vol[pidx] + gathered_vol[pidx];
          const Scalar mrel = full_vol / V_fine;

          if (mrel >= 0.5) {
            m_classifier[pidx] = PC_s;
            return;
          }

          std::vector<int> partners;

          m_particle_buckets.loop_neighbor_bucket_particles(
              bucket_idx, [&](int npidx, int) {
                if (!removed[npidx] && pidx != npidx &&
                    (m_classifier[npidx] == PC_S ||
                     m_classifier[npidx] == PC_s ||
                     m_classifier[npidx] == PC_o)) {
                  const Scalar neigh_vol = m_vol(npidx) + gathered_vol[npidx];
                  if (neigh_vol > V_fine) return false;

                  const Scalar dist =
                      (m_x.segment<3>(pidx * 3) - m_x.segment<3>(npidx * 3))
                          .norm();
                  if (dist < should_rad) {
                    partners.push_back(npidx);
                  }
                }

                return false;
              });

          if (!partners.size()) return;

          const Scalar invN = 1.0 / (Scalar)partners.size();

          const Scalar distrib_vol = full_vol * invN;
          Vec3x distrib_moment =
              (m_v.segment<3>(pidx * 3) * m_vol[pidx] + gathered_moment[pidx]) *
              invN;
          Mat3x distrib_def_grad =
              (m_Fe.block<3, 3>(pidx * 3, 0) * m_vol[pidx] +
               gathered_def_grad[pidx]) *
              invN;
          Mat3x distrib_L =
              (m_B.block<3, 3>(pidx * 3, 0) * m_vol[pidx] + gathered_L[pidx]) *
              invN;
          Mat3x distrib_b =
              (m_b.block<3, 3>(pidx * 3, 0) * m_vol[pidx] + gathered_b[pidx]) *
              invN;
          VecXx distrib_color =
              (m_components.segment(pidx * m_liquid_info.num_components,
                                    m_liquid_info.num_components) *
                   m_vol[pidx] +
               gathered_color[pidx]) *
              invN;

          for (int npidx : partners) {
            gathered_vol[npidx] += distrib_vol;
            gathered_moment[npidx] += distrib_moment;
            gathered_def_grad[npidx] += distrib_def_grad;
            gathered_L[npidx] += distrib_L;
            gathered_b[npidx] += distrib_b;
            gathered_color[npidx] += distrib_color;
          }

          removed[pidx] = true;
          m_vol(pidx) = 0.0;
          m_rest_vol(pidx) = 0.0;
          gathered_vol[pidx] = 0.0;
          gathered_moment[pidx].setZero();
          gathered_grad_p[pidx].setZero();
          gathered_def_grad[pidx].setZero();
          gathered_L[pidx].setZero();
          gathered_b[pidx].setZero();
          gathered_color[pidx].setZero();
        } else if (m_classifier[pidx] == PC_l) {
          // try downgrade level
          const Scalar full_vol = m_vol[pidx] + gathered_vol[pidx];
          if (full_vol < 1e-20) return;

          const Scalar mrel = full_vol / V_fine;

          if (mrel > 2.0) {
            m_classifier[pidx] = PC_L;
            return;
          }

          std::vector<int> partners;

          m_particle_buckets.loop_neighbor_bucket_particles(
              bucket_idx, [&](int npidx, int) {
                if (pidx != npidx && !removed[npidx] &&
                    m_classifier[npidx] == PC_s) {
                  const Scalar neigh_vol = m_vol(npidx) + gathered_vol[npidx];
                  if (neigh_vol > V_fine) return false;

                  const Scalar dist =
                      (m_x.segment<3>(pidx * 3) - m_x.segment<3>(npidx * 3))
                          .norm();
                  if (dist < should_rad) {
                    partners.push_back(npidx);
                  }
                }

                return false;
              });

          if (!partners.size()) return;

          const Scalar invN = 1.0 / (Scalar)partners.size();

          const Scalar ex_vol = std::max(full_vol - V_fine, 0.0);
          const Scalar distrib_vol = ex_vol * invN;
          const Scalar coeff = distrib_vol / full_vol;

          Vec3x distrib_moment =
              (m_v.segment<3>(pidx * 3) * m_vol[pidx] + gathered_moment[pidx]) *
              coeff;
          Mat3x distrib_def_grad =
              (m_Fe.block<3, 3>(pidx * 3, 0) * m_vol[pidx] +
               gathered_def_grad[pidx]) *
              coeff;
          Mat3x distrib_L =
              (m_B.block<3, 3>(pidx * 3, 0) * m_vol[pidx] + gathered_L[pidx]) *
              coeff;
          Mat3x distrib_b =
              (m_b.block<3, 3>(pidx * 3, 0) * m_vol[pidx] + gathered_b[pidx]) *
              coeff;
          VecXx distrib_color =
              (m_components.segment(pidx * m_liquid_info.num_components,
                                    m_liquid_info.num_components) *
                   m_vol[pidx] +
               gathered_color[pidx]) *
              coeff;

          for (int npidx : partners) {
            gathered_vol[npidx] += distrib_vol;
            gathered_moment[npidx] += distrib_moment;
            gathered_def_grad[npidx] += distrib_def_grad;
            gathered_L[npidx] += distrib_L;
            gathered_b[npidx] += distrib_b;
            gathered_color[npidx] += distrib_color;
          }

          const Scalar scaling = V_fine / full_vol;
          const Scalar rad_scaling = pow(scaling, 1.0 / 3.0);
          m_vol[pidx] *= scaling;
          m_rest_vol[pidx] *= scaling;
          m_m.segment<3>(pidx * 3) *= scaling;
          m_radius(pidx) *= rad_scaling;
          m_classifier[pidx] = PC_o;
        }
      },
      3);

  {
    LockGuard lock(m_particle_mutex);
    // gather and update
    for_each(0, num_parts, [&](int pidx) {
      if (removed[pidx] || gathered_vol[pidx] == 0.0) return;

      const Vec3x full_moment =
          m_v.segment<3>(pidx * 3) * m_vol[pidx] + gathered_moment[pidx];
      const Scalar full_vol = m_vol[pidx] + gathered_vol[pidx];

      if (full_vol == 0.0) {
        m_vol[pidx] = full_vol;
        return;
      }

      // check_isnan("unmerge J", m_J(pidx));
      // check_isnan("unmerge Fe", m_Fe.block<3, 3>(pidx * 3, 0).sum());
      // check_isnan("unmerge b trial", m_b_trial.block<3, 3>(pidx * 3,
      // 0).sum()); check_isnan("unmerge proj func", m_proj_func.segment<3>(pidx
      // * 3).sum()); check_isnan("unmerge b", m_b.block<3, 3>(pidx * 3,
      // 0).sum()); check_isnan("unmerge b plus", m_b_plus.block<3, 3>(pidx * 3,
      // 0).sum());

      const Scalar full_rad = pow(full_vol * 0.75 / M_PI, 1.0 / 3.0);
      const Scalar inv_full_vol = 1.0 / full_vol;

      // use the original J
      m_J(pidx) = m_Fe.block<3, 3>(pidx * 3, 0).determinant();

      Mat3x avg_def_grad = (m_Fe.block<3, 3>(pidx * 3, 0) * m_vol[pidx] +
                            gathered_def_grad[pidx]) *
                           inv_full_vol;

      Scalar det_avg_def_grad = avg_def_grad.determinant();
      if (det_avg_def_grad < 1e-20) {
        // deterioted case
        avg_def_grad = Mat3x::Identity();
        det_avg_def_grad = 1.0;
      }

      const Mat3x avg_L =
          (m_B.block<3, 3>(pidx * 3, 0) * m_vol[pidx] + gathered_L[pidx]) *
          inv_full_vol;
      const Mat3x avg_b =
          (m_b.block<3, 3>(pidx * 3, 0) * m_vol[pidx] + gathered_b[pidx]) *
          inv_full_vol;
      VecXx avg_c = m_components.segment(pidx * m_liquid_info.num_components,
                                         m_liquid_info.num_components) *
                        m_vol[pidx] +
                    gathered_color[pidx];
      make_gibbs_simplex(avg_c);

      m_vol[pidx] = full_vol;
      m_v.segment<3>(pidx * 3) = full_moment * inv_full_vol;
      m_m.segment<3>(pidx * 3).setConstant(full_vol * getDensity(avg_c));
      m_radius(pidx) = full_rad;
      m_Fe.block<3, 3>(pidx * 3, 0) =
          pow(m_J(pidx) / det_avg_def_grad, 1.0 / 3.0) * avg_def_grad;

      Mat3x be_bar_trial = avg_b / pow(avg_b.determinant(), 1.0 / 3.0);

      Vec3x proj_func = computeProjFunc(
          be_bar_trial, getYieldStress(avg_c), getShearModulus(avg_c),
          getViscosity(avg_c), getFlowBehaviorIndex(avg_c));
      Scalar Ie_bar = be_bar_trial.trace() / 3.0;

      m_b.block<3, 3>(pidx * 3, 0) =
          proj_func(0) * (be_bar_trial - Ie_bar * Mat3x::Identity()) +
          Ie_bar * Mat3x::Identity();
      m_b_trial.block<3, 3>(pidx * 3, 0) = be_bar_trial;
      m_Fe_plus.block<3, 3>(pidx * 3, 0) = m_Fe.block<3, 3>(pidx * 3, 0);
      m_b_plus.block<3, 3>(pidx * 3, 0) = m_b.block<3, 3>(pidx * 3, 0);
      m_proj_func.segment<3>(pidx * 3) = proj_func;
      m_B.block<3, 3>(pidx * 3, 0) = avg_L;
      m_components.segment(pidx * m_liquid_info.num_components,
                           m_liquid_info.num_components) = avg_c;

      // check_isnan("merge J", m_J(pidx));
      // check_isnan("merge Fe", m_Fe.block<3, 3>(pidx * 3, 0).sum());
      // check_isnan("merge b trial", m_b_trial.block<3, 3>(pidx * 3, 0).sum());
      // check_isnan("merge proj func", m_proj_func.segment<3>(pidx * 3).sum());
      // check_isnan("merge b", m_b.block<3, 3>(pidx * 3, 0).sum());
      // check_isnan("merge b plus", m_b_plus.block<3, 3>(pidx * 3, 0).sum());
    });

    // remove marked fluid particles
    removeEmptyParticles();

    updateOptiVolume();
  }
}

void LiquidSimulator::swapParticles(int i, int j) {
  swap<Scalar, 3>(m_x, i, j);
  swap<Scalar, 3>(m_v, i, j);
  swap<Scalar, 3>(m_m, i, j);
  std::swap(m_radius(i), m_radius(j));
  std::swap(m_J(i), m_J(j));
  std::swap(m_vol(i), m_vol(j));
  std::swap(m_rest_vol(i), m_rest_vol(j));
  std::swap(m_particle_group[i], m_particle_group[j]);
  std::swap(m_classifier[i], m_classifier[j]);
  std::swap(m_weakened[i], m_weakened[j]);

  swap<Scalar, 3>(m_B, i, j);
  swap<Scalar, 3>(m_Fe, i, j);
  swap<Scalar, 3>(m_b, i, j);
  swap<Scalar, 3>(m_b_trial, i, j);
  swap<Scalar, 3>(m_b_plus, i, j);
  swap<Scalar, 3>(m_Fe_plus, i, j);
  swap<Scalar, 3>(m_proj_func, i, j);

  swap<Scalar>(m_components, i, j, m_liquid_info.num_components);
}

void LiquidSimulator::removeEmptyParticles(int istart, bool sort) {
  const int num_parts = numParticles();

  int new_num_parts = num_parts;
  for (int i = istart; i < new_num_parts;) {
    if (m_vol(i) < 1e-20) {
      swapParticles(i, --new_num_parts);
    } else {
      ++i;
    }
  }

  if (new_num_parts < num_parts) {
    conservativeResizeParticles(new_num_parts);

    if (sort) {
      m_particle_buckets.sort(
          new_num_parts, [&](int pidx, int& i, int& j, int& k) {
            i = (int)floor((m_x(pidx * 3 + 0) - m_grid_mincorner(0)) /
                           m_bucket_size);
            j = (int)floor((m_x(pidx * 3 + 1) - m_grid_mincorner(1)) /
                           m_bucket_size);
            k = (int)floor((m_x(pidx * 3 + 2) - m_grid_mincorner(2)) /
                           m_bucket_size);
          });
    }
  }
}

bool LiquidSimulator::correctLiquidParticles() {
  const int num_fluid = numParticles();
  const Scalar dx = getDx();

  if (num_fluid == 0) return true;
  if (!m_liquid_info.correction_step) return true;

  //       std::cout << "[correctLiquidParticles 0]" << std::endl;

  m_particle_cells_single.sort(
      num_fluid, [&](int pidx, int& i, int& j, int& k) {
        Vec3x local_x = (m_x.segment<3>(pidx * 3) - m_grid_mincorner) / dx;
        i = (int)floor(local_x(0));
        j = (int)floor(local_x(1));
        k = (int)floor(local_x(2));
      });

  //        std::cout << "[correctLiquidParticles 1]" << std::endl;

  const Scalar coeff = m_liquid_info.correction_strength;

  const Scalar iD = getInverseDCoeff();

  const int correction_selector = rand() % m_liquid_info.correction_step;

  LockGuard lock(m_particle_mutex);

  m_particle_cells_single.for_each_bucket_particles_colored(
      [&](int liquid_pidx, int cell_idx) {
        if (liquid_pidx % m_liquid_info.correction_step != correction_selector)
          return;

        const Vec3x& pos = m_x.segment<3>(liquid_pidx * 3);
        Scalar dist =
            interpolateValue(pos, m_node_signed_distance, m_grid_mincorner, 0.0,
                             m_liquid_info.signed_distance_multiplier);
        if (dist > -0.75 * dx) {
          return;
        }

        const Scalar& radii = m_radius(liquid_pidx);

        Vec3x spring = Vec3x::Zero();
        m_particle_cells_single.loop_neighbor_bucket_particles(
            cell_idx, [&](int liquid_npidx, int) -> bool {
              if (liquid_pidx == liquid_npidx) return false;

              const Vec3x& np = m_x.segment<3>(liquid_npidx * 3);
              const Scalar nr = m_radius(liquid_npidx);
              const Scalar re =
                  sqrt(radii * nr) * m_liquid_info.correction_multiplier;
              const Scalar dist = (pos - np).norm();
              if (dist > re) return false;

              const Scalar w = coeff * smooth_kernel(dist * dist, re);

              if (w == 0.0) return false;

              if (dist > 1e-4 * re) {
                spring += w * (pos - np) / dist * re;
              } else {
                spring(0) += re * ScalarRand(0.0, 1.0);
                spring(1) += re * ScalarRand(0.0, 1.0);
                spring(2) += re * ScalarRand(0.0, 1.0);
              }

              return false;
            });

        Vec3x buf0 = pos + spring * m_dt;

        m_x.segment<3>(liquid_pidx * 3) = buf0;
      });

  //       std::cout << "[correctLiquidParticles 2]" << std::endl;

  return true;
}

bool LiquidSimulator::isWeakened(int pidx) const { return m_weakened[pidx]; }

bool LiquidSimulator::solidProjection() {
  const int num_parts = numParticles();
  const Scalar iD = getInverseDCoeff();
  const Scalar dx = getDx();

  // solid projection
  for_each(0, num_parts, [&](int pidx) {
    const Vec3x& pos = m_x.segment<3>(pidx * 3);
    Vec3x n = Vec3x::Zero();
    const Scalar phi = interpolateValueAndGradient(n, pos, m_node_solid_phi,
                                                   m_grid_mincorner, 3.0 * dx);

    n.normalize();

    if (m_liquid_info.solid_projection) {
      const Vec3x dpos =
          (m_v.segment<3>(pidx * 3) - get_solid_velocity(pos)) * m_dt;
      Scalar phi_now = phi + n.dot(dpos);

      if (phi < 0.0) {
        m_x.segment<3>(pidx * 3) -= phi * n;
      }
    }
  });

  return true;
}

void LiquidSimulator::terminateParticles() {
  const int num_parts = numParticles();

  for_each(0, num_parts, [&](int pidx) {
    const Vec3x& pos = m_x.segment<3>(pidx * 3);
    const Scalar phi = DistanceFieldObject::computePhi(pos, m_terminators);
    if (phi < 0.0) m_vol(pidx) = 0.0;
  });

  LockGuard lock(m_particle_mutex);

  removeEmptyParticles();
}

Scalar LiquidSimulator::getTotalVolParticles() const {
  if (!numParticles()) return 0.0;

  return m_vol.sum();
}

Scalar LiquidSimulator::getTotalVolMeshFlows() const {
  if (!m_fields.size()) return 0.0;

  VecXx field_vols(m_fields.size());
  field_vols.setZero();

  const int num_fields = (int)m_fields.size();
  for_each(0, num_fields, [&](int field_idx) {
    if (!m_fields[field_idx].has_surf_flow) return;

    std::shared_ptr<TriangularMesh> mesh =
        m_fields[field_idx].mesh_controller->getCurrentMesh();
    std::shared_ptr<TriangularMeshFlow> flow =
        m_fields[field_idx].mesh_controller->getMeshFlow();

    if (!mesh || !flow) return;

    Scalar sum_vol = 0.0;
    const int nv = mesh->nv();
    for (int i = 0; i < nv; ++i) {
      sum_vol += flow->getVertexVol(i);
    }

    field_vols[field_idx] = sum_vol;
  });

  return field_vols.sum();
}

Scalar LiquidSimulator::getTotalVolFlows() const {
  if (!m_strands.size()) return 0.0;

  VecXx strand_vols(m_strands.size());
  strand_vols.setZero();

  const int num_strands = (int)m_strands.size();
  for_each(0, num_strands, [&](int strand_idx) {
    const int num_edges = m_strands[strand_idx]->getNumEdges();
    Scalar sum_vol = 0.0;
    for (int i = 0; i < num_edges; ++i) {
      sum_vol += m_strands[strand_idx]->getFutureSurfaceFlowVolumeAtEdge(i);
    }
    strand_vols[strand_idx] = sum_vol;
  });

  return strand_vols.sum();
}

MutexType& LiquidSimulator::getParticleMutex() { return m_particle_mutex; }

MutexType& LiquidSimulator::getGridMutex() { return m_grid_mutex; }

const VecXx& LiquidSimulator::getX() const { return m_x; }
const VecXx& LiquidSimulator::getV() const { return m_v; }
const VecXx& LiquidSimulator::getM() const { return m_m; }
const VecXx& LiquidSimulator::getRadius() const { return m_radius; }
const VecXx& LiquidSimulator::getJ() const { return m_J; }
const VecXx& LiquidSimulator::getVol() const { return m_vol; }
const VecXx& LiquidSimulator::getRestVol() const { return m_rest_vol; }
const std::vector<int>& LiquidSimulator::getParticleGroup() const {
  return m_particle_group;
}
const std::vector<ParticleClassifier>& LiquidSimulator::getClassifier() const {
  return m_classifier;
}
const VecXuc& LiquidSimulator::getWeakened() const { return m_weakened; }
const VecXx& LiquidSimulator::getComponents() const { return m_components; }
const VecXx& LiquidSimulator::getProjFunc() const { return m_proj_func; }
const MatXx& LiquidSimulator::getFe() const { return m_Fe; }
const MatXx& LiquidSimulator::getb() const { return m_b; }
const MatXx& LiquidSimulator::getB() const { return m_B; }
const MatXx& LiquidSimulator::getbtrial() const { return m_b_trial; }
const MatXx& LiquidSimulator::getFePlus() const { return m_Fe_plus; }
const MatXx& LiquidSimulator::getbPlus() const { return m_b_plus; }

VecXx& LiquidSimulator::getX() { return m_x; }
VecXx& LiquidSimulator::getV() { return m_v; }
VecXx& LiquidSimulator::getM() { return m_m; }
VecXx& LiquidSimulator::getRadius() { return m_radius; }
VecXx& LiquidSimulator::getJ() { return m_J; }
VecXx& LiquidSimulator::getVol() { return m_vol; }
VecXx& LiquidSimulator::getRestVol() { return m_rest_vol; }
std::vector<int>& LiquidSimulator::getParticleGroup() {
  return m_particle_group;
}
std::vector<ParticleClassifier>& LiquidSimulator::getClassifier() {
  return m_classifier;
}
VecXuc& LiquidSimulator::getWeakened() { return m_weakened; }
VecXx& LiquidSimulator::getComponents() { return m_components; }
VecXx& LiquidSimulator::getProjFunc() { return m_proj_func; }
MatXx& LiquidSimulator::getFe() { return m_Fe; }
MatXx& LiquidSimulator::getb() { return m_b; }
MatXx& LiquidSimulator::getB() { return m_B; }
MatXx& LiquidSimulator::getbtrial() { return m_b_trial; }
MatXx& LiquidSimulator::getFePlus() { return m_Fe_plus; }
MatXx& LiquidSimulator::getbPlus() { return m_b_plus; }
};  // namespace strandsim
