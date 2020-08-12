/**
 * \copyright 2019 Yun (Raymond) Fei, 2015 Yonghao Yue
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIQUIDSIMULATOR_HH
#define LIQUIDSIMULATOR_HH

#include <vector>

#include "../Forces/Bridson/pcgsolver/sparse_matrix.hh"
#include "../Utils/Sorter.hh"
#include "../Utils/ThreadUtils.hh"
#include "DistanceFieldObject.hh"
#include "FluidScriptingController.hh"

using namespace std;

namespace strandsim {

class ElasticStrand;

enum NODE_STATE { NS_NONE, NS_FLUID, NS_SOLID };

struct LiquidInfo {
  // material-wise parameters
  VecXx liquid_density;
  VecXx surf_tension_coeff;
  VecXx liquid_bulk_modulus;
  VecXx liquid_shear_modulus;
  VecXx plastic_weaken_strain;
  VecXx plastic_yield_stress;
  VecXx flow_consistency_index;
  VecXx flow_behavior_index;
  VecXx plastic_relax_coeff;
  VecXx helmholtz_first_minima;
  VecXx helmholtz_second_minima;

  // geometric (artistic) parameters
  Scalar solid_shell_thickness;
  Scalar typical_flow_thickness;
  Scalar correction_multiplier;
  Scalar correction_strength;
  Scalar affine_stretch_damping;
  Scalar affine_rotate_damping;
  Scalar velocity_damping;
  Scalar elasto_capture_rate;
  Scalar particle_cell_multiplier;
  Scalar liquid_boundary_friction;
  Scalar theta_criterion;
  Scalar shear_pcg_criterion;
  Scalar pressure_pcg_criterion;
  Scalar levelset_young_modulus;
  Scalar mesh_flow_epsilon;
  Scalar mesh_flow_slip_length;
  Scalar mesh_flow_max_height;
  Scalar mesh_flow_critical_height;
  Scalar liquid_shear_damping;
  Scalar rest_contact_angle;
  Scalar chemical_diffusivity;
  Scalar geometric_diffusivity;
  Scalar geometric_drag_insulation;

  // numerical parameters
  int correction_step;
  int iteration_print_step;
  int surf_tension_smoothing_step;
  int pcg_max_iters;
  int num_components;
  int signed_distance_multiplier;

  // switches
  bool use_surf_tension;
  bool use_implicit_elasticity;
  bool use_implicit_pressure;
  bool solid_adhesion;
  bool solid_projection;
  bool use_liquid_capture;
  bool use_liquid_drag;
  bool use_constant_drag;
  bool use_varying_volume_fraction;
  bool solve_viscosity;
  bool solve_color_diffusion;
};

class LiquidSimulator : public FluidScriptingController,
                        public std::enable_shared_from_this<LiquidSimulator> {
 protected:
  VecXx m_x;
  VecXx m_v;
  VecXx m_m;
  VecXx m_radius;
  VecXx m_rest_vol;
  VecXx m_vol;
  VecXx m_J;
  VecXuc m_weakened;
  VecXx m_proj_func;
  VecXx m_components;

  MatXx m_Fe;
  MatXx m_b;
  MatXx m_B;
  MatXx m_b_trial;
  MatXx m_Fe_plus;
  MatXx m_b_plus;

  MatXx m_Ap;

  VecXx m_interface_points;
  vector<VecXx> m_node_signed_distance;

  vector<Mat27x2i> m_particle_nodes_x;
  vector<Mat27x2i> m_particle_nodes_y;
  vector<Mat27x2i> m_particle_nodes_z;
  vector<Mat27x2i> m_particle_nodes_p;

  vector<Mat27x4f> m_particle_weights;

  vector<Mat27x2i> m_vertex_nodes_x;
  vector<Mat27x2i> m_vertex_nodes_y;
  vector<Mat27x2i> m_vertex_nodes_z;
  vector<Mat27x2i> m_vertex_nodes_p;

  vector<Mat27x4f> m_vertex_weights;

  vector<Mat27x2i> m_inside_edge_nodes_x;
  vector<Mat27x2i> m_inside_edge_nodes_y;
  vector<Mat27x2i> m_inside_edge_nodes_z;

  vector<Mat27x3f> m_inside_edge_weights;

  typedef vector<vector<vector<pair<int, short> > > > NodeObjIndex;

  NodeObjIndex m_node_particles_x;
  NodeObjIndex m_node_particles_y;
  NodeObjIndex m_node_particles_z;
  NodeObjIndex m_node_particles_p;

  NodeObjIndex m_node_vertex_x;
  NodeObjIndex m_node_vertex_y;
  NodeObjIndex m_node_vertex_z;
  NodeObjIndex m_node_vertex_p;

  NodeObjIndex m_node_edge_x;
  NodeObjIndex m_node_edge_y;
  NodeObjIndex m_node_edge_z;

  std::vector<VecXi> m_node_pressure_neighbors;  // bucket id -> 6x2 neighbors
                                                 // of pressure nodes
  std::vector<VecXi> m_node_pp_neighbors;        // bucket id -> 18x2 pressure
                                                 // neighbors of pressure nodes

  std::vector<VecXi>
      m_node_index_pressure_x;  // bucket id -> 2x2 neighbors (left/right) to
                                // the pressure nodes
  std::vector<VecXi>
      m_node_index_pressure_y;  // bucket id -> 2x2 neighbors (down/up) to the
                                // pressure nodes
  std::vector<VecXi>
      m_node_index_pressure_z;  // bucket id -> 2x2 neighbors (back/front) to
                                // the pressure nodes

  std::vector<VecXi> m_node_index_solid_phi_x;  // bucket id -> 4x2 neighbors to
                                                // the solid phi nodes
  std::vector<VecXi> m_node_index_solid_phi_y;  // bucket id -> 4x2 neighbors to
                                                // the solid phi nodes
  std::vector<VecXi> m_node_index_solid_phi_z;  // bucket id -> 4x2 neighbors to
                                                // the solid phi nodes

  vector<int> m_particle_group;
  vector<unsigned char> m_bucket_activated;

  vector<VecXx> m_node_pos;

  vector<VecXx> m_node_solid_phi;
  vector<VecXx> m_node_liquid_phi;
  vector<VecXx> m_node_combined_phi;
  vector<VecXx> m_node_surf_tension;
  vector<VecXx> m_node_curvature_p;
  vector<VecXi> m_node_color_p;
  vector<VecXx> m_node_pressure;
  vector<VecXx> m_node_cell_solid_phi;
  vector<VecXx> m_node_components;
  vector<VecXx> m_node_elastic_vf_p;
  vector<VecXx> m_node_extra_flow_vol_p;

  vector<VecXuc> m_node_state_u;
  vector<VecXuc> m_node_state_v;
  vector<VecXuc> m_node_state_w;

  vector<VecXx> m_node_solid_vel_x;
  vector<VecXx> m_node_solid_vel_y;
  vector<VecXx> m_node_solid_vel_z;

  vector<VecXuc> m_node_liquid_valid_x;
  vector<VecXuc> m_node_liquid_valid_y;
  vector<VecXuc> m_node_liquid_valid_z;

  vector<VecXx> m_node_liquid_weight_x;  // bucket id -> node liquid weight
  vector<VecXx> m_node_liquid_weight_y;
  vector<VecXx> m_node_liquid_weight_z;

  vector<VecXx> m_node_vol_fluid_x;  // bucket id -> nodes volume
  vector<VecXx> m_node_vol_fluid_y;
  vector<VecXx> m_node_vol_fluid_z;

  vector<VecXx> m_node_vel_fluid_x;  // bucket id -> nodes velocity
  vector<VecXx> m_node_vel_fluid_y;
  vector<VecXx> m_node_vel_fluid_z;

  vector<VecXx> m_node_vel_fluid_plus_x;  // bucket id -> nodes velocity
  vector<VecXx> m_node_vel_fluid_plus_y;
  vector<VecXx> m_node_vel_fluid_plus_z;

  vector<VecXx> m_node_pressure_grad_x;  // bucket id -> nodes pressure grad
  vector<VecXx> m_node_pressure_grad_y;
  vector<VecXx> m_node_pressure_grad_z;

  //    vector< VecXx > m_node_pressure_hessian_p;

  vector<VecXx> m_node_vel_elastic_x;  // bucket id -> nodes velocity
  vector<VecXx> m_node_vel_elastic_y;
  vector<VecXx> m_node_vel_elastic_z;

  vector<VecXx> m_node_drag_coeff_x;  // bucket id -> nodes velocity
  vector<VecXx> m_node_drag_coeff_y;
  vector<VecXx> m_node_drag_coeff_z;

  vector<VecXx> m_node_vel_constraint_coeff_x;  // bucket id -> nodes velocity
  vector<VecXx> m_node_vel_constraint_coeff_y;
  vector<VecXx> m_node_vel_constraint_coeff_z;

  vector<VecXx> m_node_rhs_fluid_x;  // bucket id -> nodes forces
  vector<VecXx> m_node_rhs_fluid_y;
  vector<VecXx> m_node_rhs_fluid_z;

  vector<VecXx> m_node_mass_fluid_x;  // bucket id -> nodes mass
  vector<VecXx> m_node_mass_fluid_y;
  vector<VecXx> m_node_mass_fluid_z;

  vector<VecXx> m_node_lhs_fluid_x;  // bucket id -> nodes mass
  vector<VecXx> m_node_lhs_fluid_y;
  vector<VecXx> m_node_lhs_fluid_z;

  // assistant arrays for PCG solve
  vector<VecXx> m_node_r_x;
  vector<VecXx> m_node_r_y;
  vector<VecXx> m_node_r_z;

  vector<VecXx> m_node_z_x;
  vector<VecXx> m_node_z_y;
  vector<VecXx> m_node_z_z;

  vector<VecXx> m_node_p_x;
  vector<VecXx> m_node_p_y;
  vector<VecXx> m_node_p_z;

  vector<VecXx> m_node_q_x;
  vector<VecXx> m_node_q_y;
  vector<VecXx> m_node_q_z;

  vector<VecXx> m_node_vol_change_p;

  vector<VecXi> m_node_global_pressure_indices;
  vector<Scalar> m_pressure_rhs;
  bridson::SparseMatrix<Scalar> m_pressure_matrix;
  bridson::FixedSparseMatrix<Scalar> m_fixed_pressure_matrix;

  vector<VecXi> m_node_global_momentum_indices;
  vector<Scalar> m_momentum_rhs;
  bridson::SparseMatrix<Scalar> m_momentum_matrix;
  bridson::FixedSparseMatrix<Scalar> m_fixed_momentum_matrix;

  vector<VecXx> m_sphere_pattern;

  Vec3x m_bbx_min;
  Vec3x m_bbx_max;
  Vec3x m_grid_mincorner;

  Sorter m_particle_buckets;
  Sorter m_particle_cells;
  Sorter m_particle_cells_single;
  Sorter m_interface_buckets;
  Sorter m_vertex_buckets;
  Sorter m_edge_buckets;
  Sorter m_mesh_vertex_buckets;

  Scalar m_bucket_size;
  int m_num_nodes;
  int m_num_colers;

  LiquidInfo m_liquid_info;

  MutexType m_particle_mutex;
  MutexType m_grid_mutex;

  std::vector<DistanceFieldObject>& m_fields;
  std::vector<DistanceFieldObject>& m_sources;
  std::vector<DistanceFieldObject>& m_terminators;

  const std::vector<ElasticStrand*>& m_strands;

  vector<int> m_mesh_local_global_base;
  vector<pair<int, int> > m_mesh_global_local;

  vector<int> m_strands_local_global_base;
  vector<pair<int, int> > m_strands_global_local;
  vector<pair<int, int> > m_strands_inside_segs;
  VecXx m_strands_drag_coeffs;
  VecXx m_elastic_liquid_diffs;
  vector<unsigned char> m_strands_submerged;

  vector<Scalar> m_shooting_vol_accum;
  vector<ParticleClassifier> m_classifier;

  void sampleLiquidDistanceFields(Scalar cur_time);
  void dripFromReservoir();
  void updateParticleBoundingBox();
  void rebucketizeParticles();
  void resampleNodes();
  void updateParticleWeights();
  void updateVertexWeights();
  void updateEdgeWeights();

  void updateSolidWeights();
  void updateLiquidPhi();

  void updateOptiVolume();
  void splitLiquidParticles();
  void mergeLiquidParticles();

  void constrainLiquidVelocity();

  void terminateParticles();

  Scalar getInverseDCoeff() const;

  inline Vec3x nodePosFromBucket(int bucket_idx, int raw_node_idx,
                                 const Vec3x& offset) const;

  Vec3i getNodeHandle(int node_idx) const;

  int getNodeIndex(const Vec3i& handle) const;

  Vec3x getNodePosSolidPhi(int bucket_idx, int node_idx) const;

  Vec3x getNodePosX(int bucket_idx, int node_idx) const;

  Vec3x getNodePosY(int bucket_idx, int node_idx) const;

  Vec3x getNodePosZ(int bucket_idx, int node_idx) const;

  Vec3x getNodePosP(int bucket_idx, int node_idx) const;

  Vec3x computeProjFunc(const Mat3x& b_bar_trial,
                        const Scalar& plastic_yield_stress,
                        const Scalar& liquid_shear_modulus,
                        const Scalar& flow_consistency_index,
                        const Scalar& flow_behavior_index) const;

  void preAllocateNodes();

  void poissonSampling();

  void findNodes(const Sorter& buckets, std::vector<Mat27x2i>& particle_nodes,
                 const Vec3x& offset, const std::function<Vec3x(int)> get_pos,
                 bool activate = true);

  void expandFluidNodesMarked(int layers);

  void generateNodes();

  void connectSolidPhiNodes();

  void connectPressureNodes();

  void postAllocateNodes();

  void buildNodeParticlePairs();

  void buildNodeVertexPairs();

  void buildNodeEdgePairs();

  void transferLiquidFlowGrid();

  void solveColorDiffusion();

  template <int N>
  void buildNodePairs(const Sorter& buckets, NodeObjIndex& node_obj,
                      const vector<Eigen::Matrix<float, 27, N> >& weights,
                      const vector<Mat27x2i>& obj_node, int isel,
                      std::function<bool(int)> func = nullptr);

  void removeEmptyParticles(int istart = 0, bool sort = true);

  void swapParticles(int i, int j);

  Scalar interpolateValue(const Vec3x& pos, const std::vector<VecXx>& phi,
                          const Vec3x& phi_ori, const Scalar& default_val,
                          int subcells = 1) const;

  template <int N>
  Eigen::Matrix<Scalar, N, 1> interpolateValue(
      const Vec3x& pos, const std::vector<VecXx>& phi, const Vec3x& phi_ori,
      const Eigen::Matrix<Scalar, N, 1>& default_val, int subcells = 1) const;

  VecXx interpolateValue(const Vec3x& pos, const std::vector<VecXx>& phi,
                         const Vec3x& phi_ori, const VecXx& default_val,
                         int subcells = 1, int N = 1) const;

  Scalar interpolateValueAndGradient(Vec3x& n, const Vec3x& pos,
                                     const std::vector<VecXx>& phi,
                                     const Vec3x& phi_ori,
                                     const Scalar& default_val) const;

  Scalar getMeshFlowSlipLength() const;

  Vec3x interpolateGradient(const Vec3x& pos, const std::vector<VecXx>& phi,
                            const Vec3x& phi_ori,
                            const Scalar& default_val) const;

  template <typename T>
  void allocateNodeVectors(
      std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& node_vec_p,
      int num_subs = 1, const T& val = (T)0) const;

  template <typename T>
  void allocateNodeVectors(
      std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& node_vec_x,
      std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& node_vec_y,
      std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& node_vec_z) const;

  void applyPressureGradsFluid();

  void applyGridBasedBodyFriction(const vector<VecXx>& lhs_x,
                                  const vector<VecXx>& lhs_y,
                                  const vector<VecXx>& lhs_z);

  void solveForceExplicit(const vector<VecXx>& lhs_x,
                          const vector<VecXx>& lhs_y,
                          const vector<VecXx>& lhs_z);

  void solveForceExplicit(const vector<VecXx>& lhs_x,
                          const vector<VecXx>& lhs_y,
                          const vector<VecXx>& lhs_z,
                          const vector<VecXx>& rhs_x,
                          const vector<VecXx>& rhs_y,
                          const vector<VecXx>& rhs_z, vector<VecXx>& ret_x,
                          vector<VecXx>& ret_y, vector<VecXx>& ret_z);

  void extrapolate(vector<VecXuc>& node_valid, vector<VecXx>& node_vel);

  void computeRHS();

  void computeShearLHSDiagonal();

  void multiplyHessianMatrix(const std::vector<VecXx>& vec_x,
                             const std::vector<VecXx>& vec_y,
                             const std::vector<VecXx>& vec_z,
                             std::vector<VecXx>& ret_x,
                             std::vector<VecXx>& ret_y,
                             std::vector<VecXx>& ret_z);

  Scalar dotNodeVectors(const std::vector<VecXx>& node_vec_ax,
                        const std::vector<VecXx>& node_vec_ay,
                        const std::vector<VecXx>& node_vec_az,
                        const std::vector<VecXx>& node_vec_bx,
                        const std::vector<VecXx>& node_vec_by,
                        const std::vector<VecXx>& node_vec_bz) const;

  void particleCheckWeakened();

  Scalar solveHerschelBulkley(Scalar, Scalar, Scalar, Scalar, Scalar, Scalar);

  Scalar lengthNodeVectors(const std::vector<VecXx>& node_vec_ax,
                           const std::vector<VecXx>& node_vec_ay,
                           const std::vector<VecXx>& node_vec_az) const;

  void computeElasticLiquidDifference();

  void makeParticlesComponentConsistent();

  void makeParticlesConsistent();

  void captureBulkLiquid();

  void transferBulkLiquidMesh();

  void captureBulkLiquidMesh();

  void dripFromMesh();

  void absorbLiquidFromMesh();

 public:
  LiquidSimulator(std::vector<DistanceFieldObject>& fields_,
                  std::vector<DistanceFieldObject>& sources_,
                  std::vector<DistanceFieldObject>& terminators_,
                  const std::vector<ElasticStrand*>& strands_,
                  Scalar bucket_size_, int num_cells_, const Scalar dt_,
                  const LiquidInfo li_);

  virtual ~LiquidSimulator();

  virtual void conservativeResizeParticles(int num_particles);

  virtual void initialize();

  virtual Vec3x getNodePos(int bucket_idx, int node_idx, int r) const;
  virtual Scalar getDensityAtPos(const Vec3x& pos) const;
  virtual VecXx getComponentsAtPos(const Vec3x& pos) const;
  virtual Scalar getDensity(const VecXx& color) const;
  virtual Scalar getBulkModulus(const VecXx& color) const;
  virtual Scalar getViscosity(const VecXx& color) const;
  virtual Scalar getSolidYoungModulus() const;
  virtual Scalar getFlowBehaviorIndex(const VecXx& color) const;
  virtual Scalar getYieldStress(const VecXx& color) const;
  virtual Scalar getEpsilon(const Vec3x& pos) const;
  virtual Scalar getShearModulus(const VecXx& color) const;
  virtual Scalar getSurfTension(const VecXx& color) const;
  virtual Scalar getContactAngle() const;
  virtual int getMaxIters() const;
  virtual Scalar getCriterion() const;
  virtual Scalar getMeshFlowEpsilon() const;
  virtual Scalar getMeshFlowMaxHeight() const;
  virtual Scalar getDragInsulation() const;

  virtual bool useConstantDrag() const;

  virtual Scalar cfl();
  virtual bool map_g2p();
  virtual bool map_p2g();

  virtual bool computeSignedDistance();
  virtual bool updateSolidPhi();
  virtual bool updateDragCoeffNodes();
  virtual Vec3x addForce(const ElasticStrand*, int);
  virtual Mat3x addJacobian(const ElasticStrand*, int);
  virtual bool advectParticles();
  virtual bool addForceFluids();
  virtual bool computeWeight();
  virtual bool computePhi();
  virtual bool computeStress();
  virtual bool solvePressure();
  virtual bool solveViscosity();
  virtual bool solveDrag();
  virtual bool saveVelocity();
  virtual bool constrainVelocity();
  virtual bool prepareDataStructure();
  virtual bool processParticles();
  virtual bool solidProjection();
  virtual bool updatePlasticity();
  virtual bool applyPressure();
  virtual bool doUseDrag();
  virtual Vec3x getOrigin() const;
  virtual void getDimension(int& ni, int& nj, int& nk) const;
  virtual Scalar getDx() const;
  virtual Vec3x getParticles(int i) const;
  virtual Scalar getRadius(int i) const;
  virtual VecXx getComponents(int i) const;
  virtual bool loadParticles();
  virtual int numParticles() const;
  virtual Vec3x get_previous_velocity(const Vec3x& position) const;
  virtual Vec3x get_velocity(const Vec3x& position) const;
  virtual Scalar get_mass(const Vec3x& position) const;
  virtual Scalar get_volume(const Vec3x& position) const;
  virtual Vec3x get_future_velocity(const Vec3x& position) const;
  virtual Vec3x get_velocity(int, int) const;
  virtual Vec3x get_pressure_gradient(const Vec3x& position) const;
  virtual Mat3x get_pressure_hessian(const Vec3x& position) const;
  virtual Vec3x get_solid_gradient(const Vec3x& position) const;
  virtual Scalar get_solid_phi_gradient(const Vec3x& position, Vec3x& n) const;
  virtual Mat3x get_solid_hessian(const Vec3x& position) const;
  virtual Vec3x get_solid_velocity(const Vec3x& position) const;
  virtual Scalar get_liquid_phi(const Vec3x& position) const;
  virtual bool correctLiquidParticles();
  virtual bool isWeakened(int pidx) const;
  virtual void storeDragCoeff(int, int, const Scalar&);
  virtual bool isStrandSubmerged(int);
  virtual void findInterfacingSegments();
  virtual void updateSurfFlowConstraint();
  virtual void constrainSurfFlowVelocity();
  virtual void acceptFutureStates();
  virtual Scalar get_dense_liquid_signed_distance(const Vec3x& pos) const;
  virtual MutexType& getParticleMutex();

  virtual Scalar getTotalVolParticles() const;
  virtual Scalar getTotalVolFlows() const;
  virtual Scalar getTotalVolMeshFlows() const;
  virtual Scalar getCellSize() const;
  virtual int getNumNodes(int bucket_idx) const;
  virtual int getDefaultNumNodes() const;
  virtual MutexType& getGridMutex();
  virtual int getNumComponents() const;
  virtual Vec3i getGridDimensions() const;

  virtual const vector<VecXx>& getNodeSignedDistance() const;

  virtual const VecXx& getX() const;
  virtual const VecXx& getV() const;
  virtual const VecXx& getM() const;
  virtual const VecXx& getRadius() const;
  virtual const VecXx& getJ() const;
  virtual const VecXx& getVol() const;
  virtual const VecXx& getRestVol() const;
  virtual const std::vector<int>& getParticleGroup() const;
  virtual const std::vector<ParticleClassifier>& getClassifier() const;
  virtual const VecXuc& getWeakened() const;
  virtual const VecXx& getComponents() const;
  virtual const VecXx& getProjFunc() const;
  virtual const MatXx& getFe() const;
  virtual const MatXx& getb() const;
  virtual const MatXx& getB() const;
  virtual const MatXx& getbtrial() const;
  virtual const MatXx& getFePlus() const;
  virtual const MatXx& getbPlus() const;

  virtual VecXx& getX();
  virtual VecXx& getV();
  virtual VecXx& getM();
  virtual VecXx& getRadius();
  virtual VecXx& getJ();
  virtual VecXx& getVol();
  virtual VecXx& getRestVol();
  virtual std::vector<int>& getParticleGroup();
  virtual std::vector<ParticleClassifier>& getClassifier();
  virtual VecXuc& getWeakened();
  virtual VecXx& getComponents();
  virtual VecXx& getProjFunc();
  virtual MatXx& getFe();
  virtual MatXx& getb();
  virtual MatXx& getB();
  virtual MatXx& getbtrial();
  virtual MatXx& getFePlus();
  virtual MatXx& getbPlus();
};

}  // namespace strandsim

#endif
