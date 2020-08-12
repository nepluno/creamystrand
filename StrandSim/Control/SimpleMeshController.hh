/**
 * \copyright 2010 Breannan Smith, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SIMPLEMESHCONTROLLER_HH_
#define SIMPLEMESHCONTROLLER_HH_

#include <map>

#include "../Core/Definitions.hh"
#include "../Dynamic/MeshScriptingController.hh"
#include "ObjParser.hh"

namespace strandsim {
class TriangularMesh;
class TriangularMeshFlow;
class FluidScriptingController;
}  // namespace strandsim

namespace strandsim {
class SimpleMeshController : public MeshScriptingController {
 public:
  // static constructor
  SimpleMeshController();

  // dynamic constructor
  SimpleMeshController(const std::string& base_szfn, double startMeshTime,
                       double endMeshTime, int startMeshFrame,
                       int endMeshFrame);

  virtual ~SimpleMeshController();

  virtual void isStaticMesh(const bool i_isStaticMesh);

  virtual bool execute(bool updateLevelSet);

  const std::shared_ptr<TriangularMesh> getCurrentMesh() const {
    return m_currentMesh;
  }

  double getDefaultFrictionCoefficient() const {
    return m_defaultFrictionCoefficient;
  }

  void setDefaultFrictionCoefficient(const double frictionCoeff) {
    m_defaultFrictionCoefficient = frictionCoeff;
  }

  void getNextVelocity(Vec3xArray& v);

  virtual const std::vector<double>& getFrictionCoefficients() const {
    return m_frictionCoefficients;
  }

  virtual std::vector<double>& getFrictionCoefficients() {
    return m_frictionCoefficients;
  }

  virtual short knowsNormalSign(bool atPreviousStep, unsigned faceIndex,
                                unsigned rodIndex, unsigned vertex);

  virtual bool updateMeshNormalArea();

  bool loadNextMesh(int frame_num);

  bool updateMesh(double i_time, bool update_sub_data);

  bool initDynamicMesh(int frame_num);

  bool initMeshFlow(
      const std::shared_ptr<FluidScriptingController>& fluid_controller,
      const Scalar& init_height);

  bool initSubData();

  bool interpolateCurrentMesh(Scalar alpha);

  bool loadMesh(std::string& obj_file_name);

  bool setMesh(const TriangularMesh& trimesh);

  bool transformFlow();

  void stepFlowDynamics(const Scalar& dt);

  bool hasLevelSet() const { return false; }

  void createInitialLevelSet() {
    std::cout << "Level sets not supported by SimpleMeshController\n";
    exit(1);
  }

  LevelSet* currentLevelSet() {
    std::cout << "Level sets not supported by SimpleMeshController\n";
    exit(1);
  }

  const std::vector<bool>& getEnabledVertices() const {
    return m_enabledVertices;
  }

  virtual std::vector<bool>& getEnabledVertices() { return m_enabledVertices; }

  double getLevelSetForceThickness() const {
    std::cout << "Level sets not supported by SimpleMeshController\n";
    exit(1);
  }

  double getLevelSetForceStrength() const {
    std::cout << "Level sets not supported by SimpleMeshController\n";
    exit(1);
  }

  std::shared_ptr<TriangularMesh> getPreviousMesh() { return m_previousMesh; }

  std::shared_ptr<TriangularMesh> getNextMesh() { return m_nextMesh; }

  std::shared_ptr<TriangularMesh> getCurrentMesh() { return m_currentMesh; }

  virtual int getIFrame() { return m_currentMeshFrame; }

  virtual void setIFrame(int iframe) { m_currentMeshFrame = iframe; }

  const std::string& getBasePath() const { return m_base_szfn; }

  std::shared_ptr<TriangularMeshFlow> getMeshFlow() const { return m_meshFlow; }

 private:
  bool m_isStaticMesh;

  // Previous and next meshes are updated in this controller
  // const
  std::shared_ptr<TriangularMesh> m_previousMesh;

  std::shared_ptr<TriangularMesh> m_nextMesh;
  std::shared_ptr<TriangularMesh> m_currentMesh;

  std::shared_ptr<TriangularMeshFlow> m_meshFlow;

  std::string m_base_szfn;

  double m_startMeshTime;
  double m_endMeshTime;

  int m_currentMeshFrame;
  int m_startMeshFrame;
  int m_endMeshFrame;

  std::vector<bool> m_enabledVertices;
  std::vector<double> m_frictionCoefficients;
  double m_defaultFrictionCoefficient;
};
};  // namespace strandsim
#endif
