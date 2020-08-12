/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef __StrandSim__ProblemStepper__
#define __StrandSim__ProblemStepper__

#ifdef WIN32
#define NOMINMAX
#include <Windows.h>
#include <direct.h>
#include <shlobj_core.h>
#endif

#include <boost/generator_iterator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "../../StrandSim/Collision/CollisionUtils.hh"
#include "../../StrandSim/Control/SimpleMeshController.hh"
#include "../../StrandSim/Core/Definitions.hh"
#include "../../StrandSim/Core/ElasticStrand.hh"
#include "../../StrandSim/Dependencies/Kappas.hh"
#include "../../StrandSim/Dependencies/Twists.hh"
#include "../../StrandSim/Dynamic/DOFScriptingController.hh"
#include "../../StrandSim/Dynamic/LiquidSimulator.hh"
#include "../../StrandSim/Dynamic/StrandDynamicTraits.hh"
#include "../../StrandSim/Dynamic/StrandImplicitManager.hh"
#include "../../StrandSim/Forces/GravitationForce.hh"
#include "../../StrandSim/Forces/LevelSet.hh"
#include "../../StrandSim/Render/StrandRenderer.hh"
#include "../../StrandSim/Utils/Distances.hh"
#include "../../StrandSim/Utils/StringUtilities.hh"
#include "../../StrandSim/Utils/Timer.hh"
#include "Render/FluidsRenderer.hh"
#include "Render/TriangleMeshRenderer.hh"

using namespace strandsim;

#include <iostream>

class RodData {
 public:
  RodData(strandsim::ElasticStrand& strand, DOFScriptingController& controller,
          int num_components)
      : m_strand(strand), m_controller(controller) {
    m_renderer = new strandsim::StrandRenderer(m_strand, m_currentVertices,
                                               m_currentMaterialFrames,
                                               m_strandMutex, num_components);
    m_currentVertices.resize(m_strand.getNumVertices() * 3);
    m_currentMaterialFrames.resize((m_strand.getNumVertices() - 1) * 3);
    m_strand.dynamics().m_renderer = m_renderer;  // DEBUG
    updateDrawingData();
  }

  ~RodData() {
    if (m_renderer != NULL) {
      delete m_renderer;
      m_renderer = NULL;
    }
  }

  void update() { updateDrawingData(); }

  void updateDrawingData() {
    strandsim::LockGuard lock(m_strandMutex);

    Eigen::Map<Eigen::VectorXf> vert(&m_currentVertices[0],
                                     m_currentVertices.size());
    m_strand.getCurrentVertices(vert);

    const Vec3xArray& materialFrames1 = m_strand.getCurrentMaterialFrames1();
    int i = 0;
    for (auto mf1 = materialFrames1.begin(); mf1 != materialFrames1.end();
         ++mf1, i += 3) {
      m_currentMaterialFrames[i + 0] = (*mf1)[0];
      m_currentMaterialFrames[i + 1] = (*mf1)[1];
      m_currentMaterialFrames[i + 2] = (*mf1)[2];
    }

    m_renderer->ackGeometryChange();
  }

  void render()  // const Eigen::Matrix4f &transform )
  {
    // m_renderer->transformationMatrix() = transform;
    m_renderer->render();
  }

  Vec3x calculateObjectCenter() { return m_renderer->calculateObjectCenter(); }

  Scalar calculateObjectBoundingRadius(const Vec3d& center) {
    return m_renderer->calculateObjectBoundingRadius(center);
  }

  DOFScriptingController& getDofController() { return m_controller; }

  ElasticStrand& getStrand() { return m_strand; }

 protected:
  std::vector<float> m_currentVertices;
  std::vector<float> m_currentMaterialFrames;
  ElasticStrand& m_strand;
  DOFScriptingController& m_controller;
  strandsim::StrandRenderer* m_renderer;
  strandsim::MutexType m_strandMutex;
};

struct DumpData {
  int current_frame;
  int num_components;

  std::string fluid_fn;
  std::string rod_fn;
  std::vector<std::string> mesh_fn;

  std::vector<Vec3x> fluid_particles;
  std::vector<Scalar> fluid_radius;
  std::vector<VecXx> fluid_colors;

  std::vector<Vec4x> rod_vertices;
  std::vector<int> rod_indices;
  std::vector<Vec2x> rod_radius;
  std::vector<Scalar> rod_flow_height;
  std::vector<Scalar> rod_groups;
  std::vector<int> rod_actual;
  std::vector<VecXx> rod_colors;

  std::vector<int> rod_start_vert_indices;

  std::vector<std::vector<Vec3x> > mesh_vertices;
  std::vector<std::vector<Scalar> > mesh_flow_height;
  std::vector<std::vector<Vec3i> > mesh_indices;
  std::vector<std::vector<VecXx> > mesh_colors;
};

class ProblemStepper
    : public strandsim::StrandImplicitManager::SubStepCallback {
 public:
  virtual void executeCallback();

  explicit ProblemStepper(const std::string& name = "",
                          const std::string& desc = "");

  virtual ~ProblemStepper();

  void isSimulated(bool isSimulated) { m_isSimulated = isSimulated; }

  virtual void projectConstraint(const std::vector<ImplicitStepper*>& steppers);

  void setup(int& current_frame, int& current_check_point);

  bool step();  // take single step of problem

  void render();

  void getCenter(Vec3d& center);

  void getRadius(Scalar& radius, const Vec3d& center);

  Scalar getTime() { return m_t; }
  void setTime(Scalar t) { m_t = t; }

  Scalar getDt() { return m_dt; }
  void setDt(Scalar dt) { m_dt = dt; }

  int getSubSteps() { return m_simulation_params.m_subSteps; }
  void setSubSteps(int subSteps) { m_simulation_params.m_subSteps = subSteps; }

  // problem description :
  std::string m_problemName;
  std::string m_problemDesc;
  const std::string& ProblemName() const { return m_problemName; }
  const std::string& ProblemDescription() const { return m_problemDesc; }

  // problem options :

  virtual int LoadOptions(const char* filename) = 0;

  int LoadOptions(const std::string& filename) {
    return LoadOptions(filename.c_str());
  }

  virtual void PrintAdditionalSettings(const std::string& szfn_dir);

  void setStatLogging(int statLogging);

  void setOutputDirectory(const std::string& dir);

  // output
  void dumpMesh(DumpData* data, std::string outputdirectory, int current_frame,
                int file_width) const;
  void dumpRods(DumpData* data, std::string outputdirectory, int current_frame,
                int file_width) const;
  void dumpFluid(DumpData* data, std::string outputdirectory, int current_frame,
                 int file_width) const;

  void dumpData(std::string outputdirectory, int current_frame,
                int file_width) const;

  virtual void dumpBinaryCheckpoint(std::string outputdirectory,
                                    int current_frame, int current_checkpoint,
                                    int file_width) const;

 protected:
  int m_num_components;

  Scalar m_cfl_number;

  Scalar m_dt;

  Scalar m_t;

  bool m_isSimulated;

  virtual void setupStrands() = 0;

  virtual void setupMeshes() = 0;

  virtual void setupAfterInit(int& current_frame, int& current_check_point) = 0;

  // void AtEachTimestep();

  virtual bool executeScript(
      int total_substep_id,
      const Scalar total_substep_dt) = 0;  // execute scripting and any per-step
                                           // updating for problem

  virtual void setRodCollisionParameters(
      CollisionParameters& param,
      const ElasticStrandParameters& strand_param) = 0;

  virtual void setSimulationParameters() = 0;

  strandsim::StrandImplicitManager* m_stepper;
  std::vector<strandsim::ElasticStrand*> m_strands;

  std::vector<strandsim::ElasticStrandParameters> m_strand_parameters;

  std::vector<strandsim::CollisionParameters> m_collision_parameters;

  std::map<std::pair<int, int>, std::set<std::pair<int, int> > >
      m_collision_free;

  std::vector<RodData*> m_rodDatum;

  std::vector<std::shared_ptr<strandsim::MeshScriptingController> >
      m_meshScripting_controllers;

  std::vector<std::shared_ptr<FluidScriptingController> >
      m_fluidScripting_controllers;

  std::vector<strandsim::TriangleMeshRenderer*> m_mesh_renderers;

  std::vector<strandsim::FluidsRenderer*> m_fluids_renderers;

  std::vector<ConstraintScriptingController*> m_constraintScripting_controllers;
  SimulationParameters m_simulation_params;

  std::vector<DOFScriptingController*> m_dof_scripting_controllers;

  // hair generation
  void findOrthogonal(Vec3x& v, const Vec3x& u);
  // transform utilities
  void transformTriangleObject(TriangularMesh& triMesh,
                               const Mat3x& transformation, const Vec3x& center,
                               const Vec3x& translate);
  void transformRodRootVtx(RodData& rd, const Mat3x& transformation,
                           const Vec3x& center, const Vec3x& translate,
                           const Vec3x& scale, int vtx);

  void transformTriangleObject(TriangularMesh& triMesh,
                               const Eigen::Quaternion<double>& transformation,
                               const Vec3x& center, const Vec3x& translate);
  void transformRodRootVtx(RodData& rd,
                           const Eigen::Quaternion<double>& transformation,
                           const Vec3x& center, const Vec3x& translate,
                           const Vec3x& scale, int vtx);

  void freezeTriangleObject(TriangularMesh& triMesh);
  void translateRodVertex(RodData& rd, int vtx_id, Vec3x& translate) {
    rd.getDofController().setVertexDisplacement(
        vtx_id, translate / (Scalar)m_simulation_params.m_subSteps);
  }
};
#endif /* __StrandSim__ProblemStepper__*/
