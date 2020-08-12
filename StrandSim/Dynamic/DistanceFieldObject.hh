/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef DISTANCE_FIELD_OBJECT_HH
#define DISTANCE_FIELD_OBJECT_HH

#include "../Control/SimpleMeshController.hh"
#include "../Core/Definitions.hh"
#include "../Forces/LevelSetFwd.hh"

namespace strandsim {

enum DISTANCE_FIELD_USAGE {
  DFU_SOLID,
  DFU_SOURCE,
  DFU_TERMINATOR,

  DFU_COUNT
};

enum DISTANCE_FIELD_TYPE {
  DFT_SPHERE,
  DFT_BOX,
  DFT_CAPSULE,
  DFT_CYLINDER,
  DFT_FILE,
  DFT_SEQUENCE,

  DFT_COUNT
};

class FluidScriptingController;

struct DistanceFieldObject {
  static Scalar computePhiVel(
      const Vec3x& pos, Vec3x& vel,
      const std::vector<DistanceFieldObject>& fields,
      const std::function<bool(const DistanceFieldObject&)> selector = nullptr);

  static Scalar computePhi(
      const Vec3x& pos, const std::vector<DistanceFieldObject>& fields,
      const std::function<bool(const DistanceFieldObject&)> selector = nullptr);

  struct EMIT_INFO {
    Vec3x emit_vel;
    Scalar start;
    Scalar end;
    Scalar maxvol;
    bool enabled;
  };

  void resetDisplacement();

  DistanceFieldObject(const Vec3x& center_, const VecXx& parameter_,
                      DISTANCE_FIELD_TYPE type_, DISTANCE_FIELD_USAGE usage_,
                      const Vec3x& raxis, const double& rangle, int group_,
                      int color_index_, const double& dt, const double& dx,
                      const std::vector<DistanceFieldObject::EMIT_INFO>& emits,
                      bool inverted = false, bool clamped_by_solid = true,
                      bool has_surf_flow = false, const std::string& szfn = "",
                      const std::string& szfn_cache = "");
  virtual void advance(const double& dt, int substep_id);
  virtual void step_flow_dynamics(const double& dt);

  virtual void apply_global_rotation(const Eigen::Quaternion<double>& rot);
  virtual void apply_local_rotation(const Eigen::Quaternion<double>& rot);
  virtual void apply_translation(const Vec3x& t);
  virtual void apply_global_scaling(const Vec3x& mul);
  virtual void apply_local_scaling(const Vec3x& mul);
  virtual bool check_durations(const Scalar& cur_time, const Scalar& cur_vol,
                               Vec3x& shooting_vel) const;

  virtual void resample_internal(const std::vector<DistanceFieldObject>& solid,
                                 const Scalar& dx, const VecXx& exist,
                                 VecXx& additional, int bucket_num_nodes) const;

  void process_file_mesh(const std::string& szfn_cache);
  void init_mesh_flow(
      const std::shared_ptr<FluidScriptingController>& fluid_controller,
      const double& init_flow_height);

  bool local_bounding_box(Vec3x& bbx_low, Vec3x& bbx_high) const;

  int get_closest_face(const Vec3x& pos) const;

  Scalar compute_phi(const Vec3x& pos) const;
  Scalar compute_phi_vel(const Vec3x& pos, Vec3x& vel) const;

  DISTANCE_FIELD_TYPE type;
  DISTANCE_FIELD_USAGE usage;

  int group;

  Vec3x center;
  VecXx parameter;

  std::vector<EMIT_INFO> emits;

  std::shared_ptr<LevelSet> volume;

  Vec3x volume_origin;

  Eigen::Quaternion<double> rot;

  Vec3x future_center;
  Eigen::Quaternion<double> future_rot;
  Vec3x future_scale;

  Vec3x omega;
  Vec3x V;

  std::shared_ptr<SimpleMeshController> mesh_controller;

  Scalar sign;
  bool clamped_by_solid;
  bool has_surf_flow;

  int i_frame;
  int color_index;

  Scalar m_dx;
};

}  // namespace strandsim

#endif
