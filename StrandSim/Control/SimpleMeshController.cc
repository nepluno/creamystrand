/**
 * \copyright 2010 Breannan Smith, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "SimpleMeshController.hh"

//#include "../../Utils/Conversions.hh"

#include <iomanip>
#include <memory>

#include "../Collision/TriangularMesh.hh"
#include "../Dynamic/TriangularMeshFlow.hh"
// Duration of guessed normal sign validity, in number of frames without query
#define NORMAL_SIGN_VALIDITY 5

using namespace std;

namespace strandsim {
SimpleMeshController::SimpleMeshController(const std::string& base_szfn,
                                           double startMeshTime,
                                           double endMeshTime,
                                           int startMeshFrame, int endMeshFrame)
    : MeshScriptingController(0.,
                              (endMeshTime - startMeshTime) /
                                  (double)(endMeshFrame - startMeshFrame - 1)),
      m_base_szfn(base_szfn),
      m_isStaticMesh(false),
      m_startMeshFrame(startMeshFrame),
      m_endMeshFrame(endMeshFrame),
      m_currentMeshFrame(startMeshFrame),
      m_startMeshTime(startMeshTime),
      m_endMeshTime(endMeshTime),
      m_defaultFrictionCoefficient(0.) {
  m_currentMesh = std::make_shared<TriangularMesh>();
  m_previousMesh = std::make_shared<TriangularMesh>();
  m_nextMesh = std::make_shared<TriangularMesh>();
}

SimpleMeshController::SimpleMeshController()
    : MeshScriptingController(0., 1.),  //
      m_isStaticMesh(true),
      m_startMeshTime(0.),  //
      m_endMeshTime(0.),    //
      m_currentMeshFrame(0),
      m_startMeshFrame(0),
      m_endMeshFrame(1),
      m_defaultFrictionCoefficient(0.) {
  m_currentMesh = std::make_shared<TriangularMesh>();
  m_previousMesh = std::make_shared<TriangularMesh>();
  m_nextMesh = std::make_shared<TriangularMesh>();
}

SimpleMeshController::~SimpleMeshController() {}

void SimpleMeshController::isStaticMesh(const bool i_isStaticMesh) {
  m_isStaticMesh = i_isStaticMesh;
}

void SimpleMeshController::stepFlowDynamics(const Scalar& dt) {
  if (!m_meshFlow) return;

  m_meshFlow->step_dynamics(dt);
}

bool SimpleMeshController::transformFlow() {
  if (!m_meshFlow) return false;

  m_meshFlow->transform_base();
  return true;
}

bool SimpleMeshController::execute(bool updateLevelSet) { return true; }

bool SimpleMeshController::updateMeshNormalArea() {
  m_currentMesh->computeAdjacency(0);
  m_currentMesh->updateFaceNormalArea();
  return true;
}

bool SimpleMeshController::initSubData() {
  m_currentMesh->computeAdjacency(1);
  m_currentMesh->updateFaceNormalArea();
  m_currentMesh->computeFlattenedNeighbors();

  return true;
}

bool SimpleMeshController::initMeshFlow(
    const std::shared_ptr<FluidScriptingController>& fluid_controller,
    const Scalar& init_height) {
  m_meshFlow = std::make_shared<TriangularMeshFlow>(
      m_currentMesh, fluid_controller, init_height);
  m_currentMesh->m_associatedFlow = m_meshFlow;

  return true;
}

short SimpleMeshController::knowsNormalSign(bool atPreviousStep,
                                            unsigned faceIndex,
                                            unsigned rodIndex,
                                            unsigned vertex) {
  return 1;
}

bool SimpleMeshController::initDynamicMesh(int frame_num) {
  std::stringstream name;
  name << m_base_szfn << "." << setfill('0') << setw(5) << frame_num << ".obj";
  std::string obj_file_name = name.str().c_str();
  ObjParser objparser;
  objparser.loadTriangularMesh(obj_file_name, *m_nextMesh);
  std::cout << "Mesh Controller: Loaded init mesh " << obj_file_name << '\n';

  m_previousMesh->resizeVertices(m_nextMesh->nv());
  m_previousMesh->m_faces = m_nextMesh->m_faces;

  m_currentMesh->resizeVertices(m_nextMesh->nv());
  m_currentMesh->vertices() = m_nextMesh->vertices();
  m_currentMesh->m_faces = m_nextMesh->m_faces;

  m_previousMesh->setAssociatedController(shared_from_this());

  loadNextMesh(std::min(frame_num + 1, m_endMeshFrame - 1));
  interpolateCurrentMesh(0.0);

  return true;
}

bool SimpleMeshController::loadNextMesh(int frame_num) {
  // 1st make the previous mesh equal the next mesh as we're about to update
  // next mesh
  m_previousMesh->vertices().swap(m_nextMesh->vertices());

  // TODO : enforce that current_mesh->vertices == m_m_previousMesh->vertices()
  // here after this swap or elsewhere, i.e., in problem stepper?

  // load next mesh
  std::stringstream name;
  name << m_base_szfn << "." << setfill('0') << setw(5) << frame_num << ".obj";
  std::string obj_file_name = name.str().c_str();
  std::shared_ptr<TriangularMesh> nextFrameMesh =
      std::make_shared<TriangularMesh>();
  ObjParser objparser;
  objparser.loadTriangularMesh(obj_file_name, *nextFrameMesh);
  std::cout << "Mesh Controller: Loaded next mesh " << obj_file_name << '\n';

  for (size_t v = 0; v < m_nextMesh->nv(); ++v)
    m_nextMesh->setVertex(v, nextFrameMesh->getVertex(v));

  m_nextMesh->setAssociatedController(shared_from_this());

  return true;  // TODO: catch errors
}

bool SimpleMeshController::interpolateCurrentMesh(Scalar alpha) {
  for (size_t vtx = 0; vtx < m_previousMesh->nv(); ++vtx) {
    Vec3x currentPoint = (1. - alpha) * m_previousMesh->getVertex(vtx) +
                         alpha * m_nextMesh->getVertex(vtx);

    m_currentMesh->setDisplacement(
        vtx, currentPoint - m_currentMesh->getVertex(vtx));

    m_currentMesh->setVertex(vtx, currentPoint);
  }

  m_currentMesh->setAssociatedController(shared_from_this());

  return true;
}

void SimpleMeshController::getNextVelocity(Vec3xArray& v) {
  const Vec3xArray& x1 = m_nextMesh->vertices();
  const Vec3xArray& x0 = m_previousMesh->vertices();

  const size_t maxsize = (size_t)std::max(x0.size(), x1.size());
  v.resize(maxsize);

  for (size_t i = 0; i < maxsize; ++i) {
    v[i] = (x1[i] - x0[i]) / m_dt;
    //            std::cout << v[i] << ", " << m_dt << std::endl;
  }
}

bool SimpleMeshController::setMesh(const TriangularMesh& trimesh) {
  m_currentMesh->vertices() = trimesh.vertices();
  m_currentMesh->m_displacements = trimesh.m_displacements;
  m_currentMesh->m_faces = trimesh.m_faces;

  m_previousMesh->vertices() = trimesh.vertices();
  m_previousMesh->m_displacements = trimesh.m_displacements;
  m_previousMesh->m_faces = trimesh.m_faces;

  m_nextMesh->vertices() = trimesh.vertices();
  m_nextMesh->m_displacements = trimesh.m_displacements;
  m_nextMesh->m_faces = trimesh.m_faces;

  m_currentMesh->setAssociatedController(shared_from_this());
  m_previousMesh->setAssociatedController(shared_from_this());
  m_nextMesh->setAssociatedController(shared_from_this());

  return true;
}

bool SimpleMeshController::updateMesh(double i_time, bool update_sub_data) {
  // check if i_time is larger than next frame
  while (true) {
    if (i_time < (double)(m_currentMeshFrame - m_startMeshFrame + 1) * m_dt)
      break;

    loadNextMesh(std::min(m_currentMeshFrame + 2, m_endMeshFrame - 1));
    m_currentMeshFrame = m_currentMeshFrame + 1;
  }

  double tfrac =
      i_time / m_dt - (double)(m_currentMeshFrame - m_startMeshFrame);
  assert(tfrac >= 0. && tfrac < 1.);

  interpolateCurrentMesh(tfrac);
  m_currentMesh->updateFaceNormalArea();

  if (update_sub_data) {
    m_currentMesh->computeFlattenedNeighbors();
  }

  return true;
}

bool SimpleMeshController::loadMesh(std::string& obj_file_name) {
  // intialize all triangular meshes
  ObjParser objparser;
  objparser.loadTriangularMesh(obj_file_name, *m_currentMesh);
  // std::cout<< "Loaded mesh "<< obj_file_name << " to m_currentMesh\n";

  objparser.loadTriangularMesh(obj_file_name, *m_previousMesh);
  // std::cout<< "Loaded mesh "<< obj_file_name << " to m_previousMesh\n";

  objparser.loadTriangularMesh(obj_file_name, *m_nextMesh);
  // std::cout<< "Loaded "<< obj_file_name << " to m_nextMesh\n";

  m_currentMesh->setAssociatedController(shared_from_this());
  m_previousMesh->setAssociatedController(shared_from_this());
  m_nextMesh->setAssociatedController(shared_from_this());

  std::cout << "# Loaded mesh: \n" << obj_file_name << std::endl;

  return true;  // TODO: catch errors
}

};  // namespace strandsim
