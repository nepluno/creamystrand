/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef FLUIDSCRIPTINGCONTROLLER_HH
#define FLUIDSCRIPTINGCONTROLLER_HH

#include <vector>

#include "../Utils/ThreadUtils.hh"
#include "ScriptingController.hh"

namespace strandsim {

class ElasticStrand;

enum ParticleClassifier { PC_NON, PC_S, PC_s, PC_o, PC_l, PC_L };

class FluidScriptingController : public ScriptingController {
 public:
  struct FluidTiming {
    FluidTiming()
        : stress(0),
          datastructure(0),
          processparts(0),
          solidproj(0),
          plasticity(0),
          weight(0),
          advect(0),
          mapp2g(0),
          addforce(0),
          computephi(0),
          solvedrag(0),
          solvevisc(0),
          solvepressure(0),
          applypressure(0),
          syncflow(0),
          constraintvel(0),
          mapg2p(0) {}

    double stress;
    double datastructure;
    double processparts;
    double weight;
    double advect;
    double mapp2g;
    double addforce;
    double computephi;
    double solvedrag;
    double solvevisc;
    double syncflow;
    double solvepressure;
    double applypressure;
    double constraintvel;
    double mapg2p;
    double solidproj;
    double plasticity;

    double sum() const {
      return stress + datastructure + processparts + weight + advect + mapp2g +
             addforce + computephi + applypressure + syncflow + solvedrag +
             solvevisc + solvepressure + constraintvel + mapg2p + solidproj +
             plasticity;
    }
  };

 protected:
  std::vector<int> m_strand_offsets;
  int m_num_strand_verts;
  bool m_pressure_solved;

  std::vector<FluidTiming> m_timings;

 public:
  // get all variables
  virtual const VecXx& getX() const = 0;
  virtual const VecXx& getV() const = 0;
  virtual const VecXx& getM() const = 0;
  virtual const VecXx& getRadius() const = 0;
  virtual const VecXx& getJ() const = 0;
  virtual const VecXx& getVol() const = 0;
  virtual const VecXx& getRestVol() const = 0;
  virtual const std::vector<int>& getParticleGroup() const = 0;
  virtual const std::vector<ParticleClassifier>& getClassifier() const = 0;
  virtual const VecXuc& getWeakened() const = 0;
  virtual const VecXx& getComponents() const = 0;
  virtual const VecXx& getProjFunc() const = 0;
  virtual const MatXx& getFe() const = 0;
  virtual const MatXx& getb() const = 0;
  virtual const MatXx& getB() const = 0;
  virtual const MatXx& getbtrial() const = 0;
  virtual const MatXx& getFePlus() const = 0;
  virtual const MatXx& getbPlus() const = 0;

  virtual VecXx& getX() = 0;
  virtual VecXx& getV() = 0;
  virtual VecXx& getM() = 0;
  virtual VecXx& getRadius() = 0;
  virtual VecXx& getJ() = 0;
  virtual VecXx& getVol() = 0;
  virtual VecXx& getRestVol() = 0;
  virtual std::vector<int>& getParticleGroup() = 0;
  virtual std::vector<ParticleClassifier>& getClassifier() = 0;
  virtual VecXuc& getWeakened() = 0;
  virtual VecXx& getComponents() = 0;
  virtual VecXx& getProjFunc() = 0;
  virtual MatXx& getFe() = 0;
  virtual MatXx& getb() = 0;
  virtual MatXx& getB() = 0;
  virtual MatXx& getbtrial() = 0;
  virtual MatXx& getFePlus() = 0;
  virtual MatXx& getbPlus() = 0;

  virtual void conservativeResizeParticles(int num_particles) = 0;

  virtual void computeStrandOffsets(const std::vector<ElasticStrand*>& strands);

  FluidScriptingController(Scalar time, Scalar dt);
  virtual ~FluidScriptingController();

  virtual void initialize() = 0;

  virtual bool execute();

  virtual bool executePreStep();
  virtual bool executeMidStep();
  virtual bool executePostStep();

  virtual Scalar getShearModulus(const VecXx& color) const = 0;
  virtual Scalar getDensity(const VecXx& color) const = 0;
  virtual VecXx getComponentsAtPos(const Vec3x& pos) const = 0;
  virtual Scalar getDensityAtPos(const Vec3x& pos) const = 0;
  virtual Scalar getViscosity(const VecXx& color) const = 0;
  virtual Scalar getFlowBehaviorIndex(const VecXx& color) const = 0;
  virtual Scalar getYieldStress(const VecXx& color) const = 0;
  virtual Scalar getEpsilon(const Vec3x& pos) const = 0;
  virtual Scalar getSurfTension(const VecXx& color) const = 0;
  virtual Scalar getContactAngle() const = 0;
  virtual int getMaxIters() const = 0;
  virtual Scalar getCriterion() const = 0;

  virtual Scalar getMeshFlowMaxHeight() const = 0;
  virtual Scalar getMeshFlowSlipLength() const = 0;
  virtual Scalar getMeshFlowEpsilon() const = 0;
  virtual Vec3x computeProjFunc(const Mat3x& b_bar_trial,
                                const Scalar& plastic_yield_stress,
                                const Scalar& liquid_shear_modulus,
                                const Scalar& flow_consistency_index,
                                const Scalar& flow_behavior_index) const = 0;
  virtual Scalar cfl() = 0;
  virtual bool map_g2p() = 0;
  virtual bool map_p2g() = 0;
  virtual bool updateDragCoeffNodes() = 0;
  virtual Vec3x addForce(const ElasticStrand*, int) = 0;
  virtual Mat3x addJacobian(const ElasticStrand*, int) = 0;
  virtual bool solidProjection() = 0;
  virtual bool updatePlasticity() = 0;
  virtual bool advectParticles() = 0;
  virtual bool addForceFluids() = 0;
  virtual bool computeWeight() = 0;
  virtual bool prepareDataStructure() = 0;
  virtual bool processParticles() = 0;
  virtual bool correctLiquidParticles() = 0;
  virtual bool computePhi() = 0;
  virtual bool computeStress() = 0;
  virtual bool solvePressure() = 0;
  virtual bool applyPressure() = 0;
  virtual bool solveViscosity() = 0;
  virtual bool solveDrag() = 0;
  virtual bool constrainVelocity() = 0;
  virtual bool saveVelocity() = 0;
  virtual bool doUseDrag() = 0;
  virtual bool computeSignedDistance() = 0;
  virtual bool updateSolidPhi() = 0;
  virtual Vec3x getOrigin() const = 0;
  virtual void getDimension(int& ni, int& nj, int& nk) const = 0;
  virtual Scalar getDx() const = 0;
  virtual Vec3x getParticles(int i) const = 0;
  virtual Scalar getRadius(int i) const = 0;
  virtual bool loadParticles() = 0;
  virtual int numParticles() const = 0;
  virtual Scalar getSolidYoungModulus() const = 0;
  virtual Vec3x get_velocity(int, int) const = 0;
  virtual Vec3x get_velocity(const Vec3x& position) const = 0;
  virtual Vec3x get_future_velocity(const Vec3x& position) const = 0;
  virtual Vec3x get_previous_velocity(const Vec3x& position) const = 0;
  virtual Vec3x get_solid_velocity(const Vec3x& position) const = 0;
  virtual Scalar get_solid_phi_gradient(const Vec3x& position,
                                        Vec3x& n) const = 0;
  virtual Mat3x get_solid_hessian(const Vec3x& position) const = 0;
  virtual Vec3x get_pressure_gradient(const Vec3x& position) const = 0;
  virtual Mat3x get_pressure_hessian(const Vec3x& position) const = 0;
  virtual Scalar get_liquid_phi(const Vec3x& position) const = 0;
  virtual bool isWeakened(int pidx) const = 0;
  virtual bool useConstantDrag() const = 0;
  virtual void storeDragCoeff(int, int, const Scalar&) = 0;
  virtual bool isStrandSubmerged(int) = 0;
  virtual void findInterfacingSegments() = 0;
  virtual void updateSurfFlowConstraint() = 0;
  virtual void constrainSurfFlowVelocity() = 0;
  virtual void printTimingStats() const;
  virtual Scalar get_dense_liquid_signed_distance(const Vec3x& pos) const = 0;
  template <typename StreamT>
  void print(const FluidTiming& timings) const;

  virtual Scalar getTotalVolParticles() const = 0;
  virtual Scalar getTotalVolFlows() const = 0;
  virtual Scalar getTotalVolMeshFlows() const = 0;
  virtual VecXx getComponents(int i) const = 0;
  virtual Scalar getCellSize() const = 0;
  virtual void acceptFutureStates() = 0;
  virtual MutexType& getParticleMutex() = 0;
  virtual MutexType& getGridMutex() = 0;
  virtual int getNumNodes(int bucket_idx) const = 0;
  virtual int getDefaultNumNodes() const = 0;
  virtual int getNumComponents() const = 0;
  virtual Vec3i getGridDimensions() const = 0;
  virtual Scalar getDragInsulation() const = 0;

  virtual const std::vector<VecXx>& getNodeSignedDistance() const = 0;
  virtual Vec3x getNodePos(int bucket_idx, int node_idx, int r) const = 0;
};

FluidScriptingController::FluidTiming operator+(
    const FluidScriptingController::FluidTiming& lhs,
    const FluidScriptingController::FluidTiming& rhs);
FluidScriptingController::FluidTiming operator/(
    const FluidScriptingController::FluidTiming& lhs, const Scalar rhs);

}  // namespace strandsim

#endif
