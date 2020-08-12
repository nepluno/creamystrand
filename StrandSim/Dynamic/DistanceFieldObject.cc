/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "DistanceFieldObject.hh"

#include <iomanip>

#include "../Control/Capsule.hh"
#include "../Control/Icosphere.hh"
#include "../Control/ObjParser.hh"
#include "../Control/RoundCornerBox.hh"
#include "../Control/RoundCylinder.hh"
#include "../Forces/Bridson/array3.hh"
#include "../Forces/Bridson/array3_utils.hh"
#include "../Forces/LevelSet.hh"
#include "../Utils/MathUtilities.hh"
#include "../Utils/Sorter.hh"
#include "../Utils/ThreadUtils.hh"

namespace strandsim {
using namespace bridson;

inline Scalar sphere_phi(const Vec3x& position, const Vec3x& centre,
                         Scalar radius) {
  return ((position - centre).norm() - radius);
}

inline void sphere_phi_bbx(const Vec3x& centre, Scalar radius, Vec3x& bbx_low,
                           Vec3x& bbx_high) {
  bbx_low = centre - Vec3x(radius, radius, radius);
  bbx_high = centre + Vec3x(radius, radius, radius);
}

inline Scalar capsule_phi(const Vec3x& position, const Vec3x& centre,
                          const Scalar& radius, const Scalar& halflength) {
  Vec3x a = centre - Vec3x(halflength, 0, 0);
  Vec3x pa = position - a;
  Vec3x ba = Vec3x(2.0 * halflength, 0, 0);
  Scalar h = clamp(pa.dot(ba) / (4.0 * halflength * halflength), 0.0, 1.0);
  return (pa - ba * h).norm() - radius;
}

inline void capsule_phi_bbx(const Vec3x& centre,
                            const Eigen::Quaternion<Scalar>& rot,
                            const Scalar& radius, const Scalar& halflength,
                            Vec3x& bbx_low, Vec3x& bbx_high) {
  const Vec3x p0 = rot * Vec3x(-halflength, 0, 0) + centre;
  const Vec3x p1 = rot * Vec3x(halflength, 0, 0) + centre;

  bbx_low = Vec3x(std::min(p0(0), p1(0)), std::min(p0(1), p1(1)),
                  std::min(p0(2), p1(2))) -
            Vec3x(radius, radius, radius);
  bbx_high = Vec3x(std::max(p0(0), p1(0)), std::max(p0(1), p1(1)),
                   std::max(p0(2), p1(2))) +
             Vec3x(radius, radius, radius);
}

inline Scalar box_phi(const Vec3x& position, const Vec3x& centre,
                      const Vec3x& expand, const Scalar& radius) {
  Scalar dx = fabs(position[0] - centre[0]) - expand[0];
  Scalar dy = fabs(position[1] - centre[1]) - expand[1];
  Scalar dz = fabs(position[2] - centre[2]) - expand[2];
  Scalar dax = std::max(dx, 0.0);
  Scalar day = std::max(dy, 0.0);
  Scalar daz = std::max(dz, 0.0);
  return std::min(std::max(std::max(dx, dy), dz), 0.0) +
         sqrt(dax * dax + day * day + daz * daz) - radius;
}

inline void box_phi_bbx(const Vec3x& centre,
                        const Eigen::Quaternion<Scalar>& rot,
                        const Vec3x& expand, const Scalar& radius,
                        Vec3x& bbx_low, Vec3x& bbx_high) {
  bbx_low = bbx_high = centre;

  for (int r = 0; r < 2; ++r)
    for (int s = 0; s < 2; ++s)
      for (int t = 0; t < 2; ++t) {
        Vec3x p =
            rot * Vec3x(r ? expand(0) : -expand(0), s ? expand(1) : -expand(1),
                        t ? expand(2) : -expand(2)) +
            centre;
        bbx_low = Vec3x(std::min(bbx_low(0), p(0)), std::min(bbx_low(1), p(1)),
                        std::min(bbx_low(2), p(2)));
        bbx_high =
            Vec3x(std::max(bbx_high(0), p(0)), std::max(bbx_high(1), p(1)),
                  std::max(bbx_high(2), p(2)));
      }

  bbx_low -= Vec3x(radius, radius, radius);
  bbx_high += Vec3x(radius, radius, radius);
}

inline Scalar cylinder_phi(const Vec3x& position, const Vec3x& centre,
                           const Scalar& radius_ext, const Scalar& radius_cor,
                           const Scalar& h) {
  Vec3x p = position - centre;
  Vec2x d =
      Vec2x(Vec2x(p(0), p(2)).norm() - radius_ext + radius_cor, fabs(p(1)) - h);
  return std::min(std::max(d(0), d(1)), 0.0) +
         Vec2x(std::max(d(0), 0.0), std::max(d(1), 0.0)).norm() - radius_cor;
}

inline void cylinder_phi_bbx(const Vec3x& centre,
                             const Eigen::Quaternion<Scalar>& rot,
                             const Scalar& radius_ext, const Scalar& radius_cor,
                             const Scalar& h, Vec3x& bbx_low, Vec3x& bbx_high) {
  Vec3x ptb = rot * Vec3x(0, h + radius_cor, 0) + centre;
  Vec3x pta = rot * Vec3x(0, -h - radius_cor, 0) + centre;
  Vec3x a = ptb - pta;
  Scalar da = a.dot(a);

  Vec3x db = (radius_ext + radius_cor) * Vec3x(sqrt(1.0 - a(0) * a(0) / da),
                                               sqrt(1.0 - a(1) * a(1) / da),
                                               sqrt(1.0 - a(2) * a(2) / da));
  bbx_low = Vec3x(std::min(pta(0), ptb(0)), std::min(pta(1), ptb(1)),
                  std::min(pta(2), ptb(2))) -
            db;
  bbx_high = Vec3x(std::max(pta(0), ptb(0)), std::max(pta(1), ptb(1)),
                   std::max(pta(2), ptb(2))) +
             db;

  Vec3x expansion = (bbx_high - bbx_low) * 0.5;
  bbx_low = centre - expansion;
  bbx_high = centre + expansion;
}

void DistanceFieldObject::resetDisplacement() {
  if (!mesh_controller) return;

  mesh_controller->getCurrentMesh()->zeroDisplacement();
  mesh_controller->getNextMesh()->zeroDisplacement();
  mesh_controller->getPreviousMesh()->zeroDisplacement();
}

void DistanceFieldObject::apply_global_scaling(const Vec3x& mul) {
  future_scale.array() *= mul.array();
  future_center.array() *= mul.array();
}

void DistanceFieldObject::apply_local_scaling(const Vec3x& mul) {
  future_scale.array() *= mul.array();
}

DistanceFieldObject::DistanceFieldObject(
    const Vec3x& center_, const VecXx& parameter_, DISTANCE_FIELD_TYPE type_,
    DISTANCE_FIELD_USAGE usage_, const Vec3x& raxis, const double& rangle,
    int group_, int color_index_, const double& dt, const double& dx,
    const std::vector<DistanceFieldObject::EMIT_INFO>& emits_, bool inverted,
    bool clamped_by_solid_, bool has_surf_flow_, const std::string& szfn,
    const std::string& szfn_cache)
    : type(type_),
      usage(usage_),
      group(group_),
      color_index(color_index_),
      m_dx(dx),
      center(center_),
      parameter(parameter_),
      emits(emits_),
      rot(Eigen::AngleAxis<double>(rangle, raxis)),
      sign(inverted ? -1.0 : 1.0),
      clamped_by_solid(clamped_by_solid_),
      has_surf_flow(has_surf_flow_),
      i_frame(0),
      future_scale(Vec3x::Ones()) {
  V.setZero();
  omega.setZero();

  future_center = center;
  future_rot = rot;

  if (type_ == DFT_SEQUENCE) {
    mesh_controller = std::make_shared<SimpleMeshController>(
        szfn, parameter(0), parameter(1), (int)parameter(2), (int)parameter(3));
    mesh_controller->initDynamicMesh((int)parameter(2));

    std::stringstream oss;
    oss << mesh_controller->getBasePath() << "." << dt << "." << m_dx << "."
        << std::setfill('0') << std::setw(5) << i_frame << ".cache";
    process_file_mesh(oss.str());

    mesh_controller->getCurrentMesh()->buildTriangleIndexer(m_dx * 2.0);
  } else {
    mesh_controller = std::make_shared<SimpleMeshController>();

    TriangularMesh newMesh;

    newMesh.setAssociatedController(mesh_controller);

    switch (type_) {
      case DFT_BOX: {
        RoundCornerBox(32, Vec3x(parameter(0), parameter(1), parameter(2)),
                       parameter(3), &newMesh, inverted);
      } break;
      case DFT_SPHERE: {
        Icosphere(4, parameter(0), &newMesh, inverted);
      } break;
      case DFT_CAPSULE: {
        Capsule(128, parameter(0), parameter(1), &newMesh, inverted);

        break;
      }
      case DFT_CYLINDER: {
        RoundCylinder(32, 8, parameter(0), parameter(1), parameter(2), &newMesh,
                      inverted);

        break;
      }
      case DFT_FILE: {
        ObjParser(szfn, &newMesh, inverted);
        newMesh.scale(Vec3x::Constant(parameter(0)));

        break;
      }
      default:
        break;
    }

    mesh_controller->setMesh(newMesh);

    if (type_ == DFT_FILE) {
      process_file_mesh(szfn_cache);
    }
  }

  if (usage_ == DFU_SOLID) {
    if (has_surf_flow)
      mesh_controller->initSubData();
    else
      mesh_controller->updateMeshNormalArea();
  }
}

void DistanceFieldObject::init_mesh_flow(
    const std::shared_ptr<FluidScriptingController>& fluid_controller,
    const double& init_flow_height) {
  if (usage == DFU_SOLID && has_surf_flow) {
    mesh_controller->initMeshFlow(fluid_controller, init_flow_height);
  }
}

int DistanceFieldObject::get_closest_face(const Vec3x& pos) const {
  if (type == DFT_SEQUENCE) {
    return mesh_controller->getCurrentMesh()->getClosestTriangle(pos);
  } else {
    return -1;
  }
}

void DistanceFieldObject::process_file_mesh(const std::string& szfn_cache) {
  std::shared_ptr<TriangularMesh> mesh = mesh_controller->getCurrentMesh();

  std::vector<unsigned int> tri_indices(mesh->nf());
  for (unsigned int i = 0; i < mesh->nf(); ++i) {
    tri_indices[i] = i;
  }

  // pre-process mesh

  Vec3x bbx_min, bbx_max;
  mesh->computeBBox();
  bbx_min = mesh->getBBoxCenter() - mesh->getBBoxExtent();
  bbx_max = mesh->getBBoxCenter() + mesh->getBBoxExtent();

  bbx_min -= Vec3x::Constant(m_dx * 3.0);
  bbx_max += Vec3x::Constant(m_dx * 3.0);

  volume_origin = bbx_min;

  volume = std::make_shared<LevelSet>();

  Vec3x extend = (bbx_max - bbx_min) / m_dx;

  int nx = (int)std::ceil(extend(0));
  int ny = (int)std::ceil(extend(1));
  int nz = (int)std::ceil(extend(2));

  Vec3xArray v;
  mesh_controller->getNextVelocity(v);

  // check if cache exist
  if (!szfn_cache.empty()) {
    std::ifstream ifs(szfn_cache, std::ios::binary);

    if (ifs.good()) {
      // check file size
      size_t i_fz_ask = volume->getFileSize(
          volume_origin, volume_origin + Vec3x(nx, ny, nz) * m_dx * 0.5, m_dx);

      std::ifstream ifs_test(szfn_cache,
                             std::ifstream::ate | std::ifstream::binary);
      size_t i_fz_get = ifs_test.tellg();

      if (i_fz_ask == i_fz_get) {
        // read from file directly
        std::cout << "[read from cache " << szfn_cache << "]" << std::endl;

        volume->buildLevelSetFromFile(ifs, mesh->getFaces(), tri_indices,
                                      mesh->vertices(), v, volume_origin, m_dx,
                                      nx, ny, nz);
        ifs.close();
        return;
      } else {
        std::cout << "[cache " << szfn_cache << " has invalid size ("
                  << i_fz_get << ", " << i_fz_ask << ")]" << std::endl;
      }
    }
  }

  if (sign < 0.0) mesh->invert();

  volume->buildLevelSet(mesh->getFaces(), tri_indices, mesh->vertices(), v,
                        volume_origin, m_dx, nx, ny, nz);

  if (sign < 0.0) mesh->invert();

  if (!szfn_cache.empty()) {
    volume->writeFileAsync(szfn_cache);
  }
}

void DistanceFieldObject::apply_global_rotation(
    const Eigen::Quaternion<double>& rot) {
  future_rot = rot * future_rot;
}

void DistanceFieldObject::apply_local_rotation(
    const Eigen::Quaternion<double>& rot) {
  future_rot = future_rot * rot;
}

void DistanceFieldObject::apply_translation(const Vec3x& t) {
  future_center += t;
}

void DistanceFieldObject::step_flow_dynamics(const double& dt) {
  if (usage == DFU_SOLID && has_surf_flow) {
    mesh_controller->stepFlowDynamics(dt);
  }
}

void DistanceFieldObject::advance(const double& dt, int substep_id) {
  mesh_controller->getCurrentMesh()->saveDisplacement();

  V = (future_center - center) / dt;
  Eigen::Quaternion<double> q = future_rot * rot.conjugate();
  double len = q.vec().norm();
  if (len > 0.0) {
    double angle = 2.0 * atan2(len, q.w());
    omega = q.vec() / len * angle / dt;
  } else {
    omega.setZero();
  }

  center = future_center;
  rot = future_rot;

  parameter(0) *= future_scale(0);
  parameter(1) *= future_scale(1);
  parameter(2) *= future_scale(2);

  i_frame++;

  mesh_controller->getCurrentMesh()->m_stored_vertices =
      mesh_controller->getCurrentMesh()->m_vertices;

  LockGuard lock(mesh_controller->getCurrentMesh()->getGeometryMutex());

  if (type == DFT_SEQUENCE) {
    mesh_controller->updateMesh(i_frame * dt, has_surf_flow);
    mesh_controller->getCurrentMesh()->transform(future_rot, Vec3x::Zero(),
                                                 future_center, future_scale);
    std::stringstream oss;
    oss << mesh_controller->getBasePath() << "." << dt << "." << m_dx << "."
        << std::setfill('0') << std::setw(5) << i_frame << "."
        << std::setfill('0') << std::setw(5) << substep_id << ".cache";
    process_file_mesh(oss.str());
  } else {
    mesh_controller->getCurrentMesh()->m_vertices =
        mesh_controller->getPreviousMesh()->m_vertices;
    mesh_controller->getCurrentMesh()->transform(future_rot, Vec3x::Zero(),
                                                 future_center, future_scale);
  }

  if (usage == DFU_SOLID && has_surf_flow) {
    mesh_controller->getCurrentMesh()->updateFaceNormalArea();
    mesh_controller->transformFlow();
  }

  //        future_scale.setOnes();
}

bool DistanceFieldObject::check_durations(const Scalar& cur_time,
                                          const Scalar& cur_vol,
                                          Vec3x& shooting_vel) const {
  bool ret = false;
  for (auto& dur : emits) {
    ret =
        (cur_time >= dur.start && cur_time <= dur.end) && cur_vol < dur.maxvol;
    if (ret) {
      shooting_vel = dur.emit_vel;
      break;
    }
  }

  return ret;
}

Scalar DistanceFieldObject::compute_phi_vel(const Vec3x& pos,
                                            Vec3x& vel) const {
  Scalar phi;

  if (type == DFT_SEQUENCE) {
    Vec3x dx = pos - center;
    Eigen::Quaternion<Scalar> p0(0.0, dx(0), dx(1), dx(2));
    Eigen::Quaternion<Scalar> irot = future_rot.conjugate();
    Vec3x rotp = (irot * p0 * irot.inverse()).vec();
    phi = volume->getLevelSetValueVelocity(rotp, vel);
    //#pragma omp critical
    //            {
    //                std::cout << vel << std::endl;
    //            }
    //            vel += V + omega.cross(dx);
  } else {
    phi = compute_phi(pos);

    Vec3x dx = pos - center;

    vel = V + omega.cross(dx);
  }

  return phi;
}

Scalar DistanceFieldObject::compute_phi(const Vec3x& pos) const {
  Scalar phi = 0.0;

  Vec3x dx = pos - center;

  switch (type) {
    case DFT_BOX: {
      Eigen::Quaternion<Scalar> p0(0.0, dx(0), dx(1), dx(2));
      Eigen::Quaternion<Scalar> irot = future_rot.conjugate();
      Vec3x rotp = (irot * p0 * irot.inverse()).vec() + center;
      phi = sign * box_phi(rotp, center,
                           Vec3x(parameter(0), parameter(1), parameter(2)),
                           parameter(3));
      break;
    }

    case DFT_SPHERE: {
      phi = sign * sphere_phi(pos, center, parameter(0));
      break;
    }

    case DFT_CAPSULE: {
      Eigen::Quaternion<Scalar> p0(0.0, dx(0), dx(1), dx(2));
      Eigen::Quaternion<Scalar> irot = future_rot.conjugate();
      Vec3x rotp = (irot * p0 * irot.inverse()).vec() + center;
      phi = sign * capsule_phi(rotp, center, parameter(0), parameter(1));
      break;
    }

    case DFT_CYLINDER: {
      Eigen::Quaternion<Scalar> p0(0.0, dx(0), dx(1), dx(2));
      Eigen::Quaternion<Scalar> irot = future_rot.conjugate();
      Vec3x rotp = (irot * p0 * irot.inverse()).vec() + center;
      phi = sign * cylinder_phi(rotp, center, parameter(0), parameter(1),
                                parameter(2));
      break;
    }

    case DFT_FILE:
    case DFT_SEQUENCE: {
      Eigen::Quaternion<Scalar> p0(0.0, dx(0), dx(1), dx(2));
      Eigen::Quaternion<Scalar> irot = future_rot.conjugate();
      Vec3x rotp = (irot * p0 * irot.inverse()).vec();
      rotp.array() /= future_scale.array();
      phi = sign * volume->getLevelSetValue(rotp);
      break;
    }
    default:
      break;
  }

  return phi;
}

bool DistanceFieldObject::local_bounding_box(Vec3x& bbx_low,
                                             Vec3x& bbx_high) const {
  switch (type) {
    case DFT_CAPSULE:
      capsule_phi_bbx(center, rot, parameter(0), parameter(1), bbx_low,
                      bbx_high);
      break;
    case DFT_SPHERE:
      sphere_phi_bbx(center, parameter(0), bbx_low, bbx_high);
      break;
    case DFT_BOX:
      box_phi_bbx(center, rot, Vec3x(parameter(0), parameter(1), parameter(2)),
                  parameter(3), bbx_low, bbx_high);
      break;
    case DFT_CYLINDER:
      cylinder_phi_bbx(center, rot, parameter(0), parameter(1), parameter(2),
                       bbx_low, bbx_high);
      break;
    case DFT_FILE:
    case DFT_SEQUENCE:
      mesh_controller->getCurrentMesh()->computeBBox();
      bbx_low = mesh_controller->getCurrentMesh()->m_bbox_center -
                mesh_controller->getCurrentMesh()->m_bbox_extent;
      bbx_high = mesh_controller->getCurrentMesh()->m_bbox_center +
                 mesh_controller->getCurrentMesh()->m_bbox_extent;
      break;
    default:
      break;
  }

  return sign < 0.0;
}

void DistanceFieldObject::resample_internal(
    const std::vector<DistanceFieldObject>& solid, const Scalar& dx,
    const VecXx& exist, VecXx& additional, int bucket_num_nodes) const {
  Vec3x bbx_low = Vec3x::Zero();
  Vec3x bbx_high = Vec3x::Zero();

  local_bounding_box(bbx_low, bbx_high);

  bbx_low -= Vec3x(dx, dx, dx);
  bbx_high += Vec3x(dx, dx, dx);

  const Scalar bucket_size = dx * bucket_num_nodes;

  Vec3x dxyz = (bbx_high - bbx_low) / bucket_size;
  if (dxyz(0) <= 0.0 || dxyz(1) <= 0.0 || dxyz(2) <= 0.0) return;

  Vec3i nxyz = Vec3i((int)std::ceil(dxyz(0)), (int)std::ceil(dxyz(1)),
                     (int)std::ceil(dxyz(2)));

  // build buckets
  Sorter sorter(nxyz(0), nxyz(1), nxyz(2));

  const int num_parts = (int)exist.size() / 3;

  sorter.sort(num_parts, [&](int pidx, int& i, int& j, int& k) {
    i = (int)std::floor((exist(pidx * 3 + 0) - bbx_low(0)) / bucket_size);
    j = (int)std::floor((exist(pidx * 3 + 1) - bbx_low(1)) / bucket_size);
    k = (int)std::floor((exist(pidx * 3 + 2) - bbx_low(2)) / bucket_size);
  });

  std::vector<std::vector<Vec3x> > bucket_samples(nxyz(0) * nxyz(1) * nxyz(2));

  sorter.for_each_bucket([&](int bucket_idx) {
    Vec3i handle = sorter.bucket_handle(bucket_idx);
    Vec3x base_pos =
        bbx_low + Vec3x(handle(0), handle(1), handle(2)) * bucket_size;

    std::vector<unsigned char> counter(
        bucket_num_nodes * bucket_num_nodes * bucket_num_nodes, 0U);

    sorter.get_bucket(bucket_idx, [&](int pidx) {
      int i = (int)std::floor((exist(pidx * 3 + 0) - base_pos(0)) / dx);
      int j = (int)std::floor((exist(pidx * 3 + 1) - base_pos(1)) / dx);
      int k = (int)std::floor((exist(pidx * 3 + 2) - base_pos(2)) / dx);

      if (i >= 0 && i < bucket_num_nodes && j >= 0 && j < bucket_num_nodes &&
          k >= 0 && k < bucket_num_nodes) {
        int local_cell_idx =
            k * bucket_num_nodes * bucket_num_nodes + j * bucket_num_nodes + i;
        counter[local_cell_idx] = 1U;
      }
    });

    for (int k = 0; k < bucket_num_nodes; ++k)
      for (int j = 0; j < bucket_num_nodes; ++j)
        for (int i = 0; i < bucket_num_nodes; ++i) {
          int local_cell_idx = k * bucket_num_nodes * bucket_num_nodes +
                               j * bucket_num_nodes + i;

          // ignore occupied cells
          if (counter[local_cell_idx]) continue;

          // sample if empty
          Vec3x cell_min_corner = base_pos + Vec3x(i, j, k) * dx;
          Vec3x cell_max_corner = cell_min_corner + Vec3x(dx, dx, dx);

          Vec3x p = Vec3x(ScalarRand(cell_min_corner(0), cell_max_corner(0)),
                          ScalarRand(cell_min_corner(1), cell_max_corner(1)),
                          ScalarRand(cell_min_corner(2), cell_max_corner(2)));

          // test against source levelset
          Scalar phi = compute_phi(p);
          if (phi > 0.0) continue;

          // test against solid levelset
          if (clamped_by_solid) {
            bool inside_solid = false;
            for (const DistanceFieldObject& s : solid) {
              if (s.compute_phi(p) < 0.0) {
                inside_solid = true;
                break;
              }
            }

            if (inside_solid) continue;
          }

          // accept new sample
          bucket_samples[bucket_idx].push_back(p);
        }
  });

  parallel_concatenate(additional, bucket_samples);
}

Scalar DistanceFieldObject::computePhiVel(
    const Vec3x& pos, Vec3x& vel,
    const std::vector<DistanceFieldObject>& fields,
    const std::function<bool(const DistanceFieldObject&)> selector) {
  Scalar min_phi = std::numeric_limits<Scalar>::max();
  Vec3x min_vel = Vec3x::Zero();
  for (auto dfptr : fields) {
    if (selector && !selector(dfptr)) continue;

    Vec3x v;
    Scalar phi = dfptr.compute_phi_vel(pos, v);
    if (phi < min_phi) {
      min_phi = phi;
      min_vel = v;
    }
  }

  vel = min_vel;

  return min_phi;
}

Scalar DistanceFieldObject::computePhi(
    const Vec3x& pos, const std::vector<DistanceFieldObject>& fields,
    const std::function<bool(const DistanceFieldObject&)> selector) {
  Scalar min_phi = std::numeric_limits<Scalar>::max();
  for (auto dfptr : fields) {
    if (selector && !selector(dfptr)) continue;

    Scalar phi = dfptr.compute_phi(pos);
    if (phi < min_phi) {
      min_phi = phi;
    }
  }

  return min_phi;
}

}  // namespace strandsim
