/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "FluidScriptingController.hh"

#include "../Core/ElasticStrand.hh"
#include "../Forces/FluidDragForce.hh"
#include "../Utils/TextLog.hh"
#include "../Utils/TimeUtils.hh"
#include "../Utils/Timer.hh"

namespace strandsim {

FluidScriptingController::FluidScriptingController(double time, double dt)
    : ScriptingController(time, dt) {}

FluidScriptingController::~FluidScriptingController() {
  FluidDragForce::setScriptingController(NULL);
}

bool FluidScriptingController::execute() {
  bool a = executePreStep();
  bool b = executePostStep();
  return a & b;
}

bool FluidScriptingController::executePreStep() {
  Timer timer("fluid", false);
  FluidTiming timing;
  bool res = true;

  std::cout << "[fluid: process particles]" << std::endl;
  timer.restart();
  res &= processParticles();
  timing.processparts = timer.elapsed();

  std::cout << "[fluid: prepare data structure]" << std::endl;
  timer.restart();
  res &= prepareDataStructure();
  timing.datastructure = timer.elapsed();

  std::cout << "[fluid: update solid phi]" << std::endl;
  timer.restart();
  res &= updateSolidPhi();

  std::cout << "[fluid: compute signed distance]" << std::endl;
  res &= computeSignedDistance();

  std::cout << "[fluid: correct liquid particles]" << std::endl;
  res &= correctLiquidParticles();

  std::cout << "[fluid: compute weight]" << std::endl;
  res &= computeWeight();
  timing.weight = timer.elapsed();

  std::cout << "[fluid: compute phi]" << std::endl;
  timer.restart();
  res &= computePhi();
  timing.computephi = timer.elapsed();

  std::cout << "[fluid: map_p2g]" << std::endl;
  timer.restart();
  res &= map_p2g();
  timing.mapp2g = timer.elapsed();

  m_timings.push_back(timing);
  return res;
}

bool FluidScriptingController::executeMidStep() {
  bool res = true;
  Timer timer("fluid", false);
  FluidTiming timing = m_timings.back();
  m_timings.pop_back();

  res &= updateDragCoeffNodes();

  std::cout << "[fluid: add force]" << std::endl;
  timer.restart();
  res &= addForceFluids();
  timing.addforce = timer.elapsed();

  std::cout << "[fluid: solve viscosity before pressure]" << std::endl;
  timer.restart();
  res &= solveViscosity();
  timing.solvevisc = timer.elapsed();

  std::cout << "[fluid: solve pressure]" << std::endl;
  timer.restart();
  res &= solvePressure();
  timing.solvepressure = timer.elapsed();

  std::cout << "[fluid: apply new pressure]" << std::endl;
  timer.restart();
  res &= applyPressure();
  timing.applypressure = timer.elapsed();

  std::cout << "[fluid: solve viscosity after pressure]" << std::endl;
  timer.restart();
  res &= solveViscosity();
  timing.solvevisc = timer.elapsed();

  std::cout << "[fluid: constrain depth-averaged flow]" << std::endl;
  timer.restart();
  constrainSurfFlowVelocity();
  timing.syncflow += timer.elapsed();

  std::cout << "[fluid: constrain velocity]" << std::endl;
  timer.restart();
  res &= constrainVelocity();
  timing.constraintvel = timer.elapsed();

  std::cout << "[fluid: map_g2p]" << std::endl;
  timer.restart();
  res &= map_g2p();
  timing.mapg2p = timer.elapsed();

  acceptFutureStates();

  m_timings.push_back(timing);
  return res;
}

bool FluidScriptingController::executePostStep() {
  bool res = true;

  Timer timer("fluid", false);
  FluidTiming timing = m_timings.back();
  m_timings.pop_back();

  std::cout << "[fluid: advect particles]" << std::endl;
  timer.restart();
  res &= advectParticles();
  timing.advect = timer.elapsed();

  std::cout << "[fluid: solid projection]" << std::endl;
  timer.restart();
  res &= solidProjection();
  timing.solidproj = timer.elapsed();

  std::cout << "[fluid: update plasticity]" << std::endl;
  timer.restart();
  res &= updatePlasticity();
  timing.plasticity = timer.elapsed();

  m_timings.push_back(timing);

  return res;
}

void FluidScriptingController::computeStrandOffsets(
    const std::vector<ElasticStrand*>& strands) {
  m_strand_offsets.resize(strands.size());

  int k = 0;
  for (int i = 0; i < strands.size(); ++i) {
    m_strand_offsets[i] = k;
    k += strands[i]->getNumVertices();
  }
  m_num_strand_verts = k;
}

void FluidScriptingController::printTimingStats() const {
  FluidTiming tot;

  for (auto stepT = m_timings.begin(); stepT != m_timings.end(); ++stepT) {
    tot = tot + *stepT;
  }

  print<InfoStream>(tot / (Scalar)m_timings.size());

  const Scalar vol_parts = getTotalVolParticles();
  const Scalar vol_flows = getTotalVolFlows();
  const Scalar vol_mesh_flows = getTotalVolMeshFlows();

  const int num_parts = numParticles();
  const Vec3i dim_grid = getGridDimensions();

  std::cout << "[fluid # particles: " << num_parts << ", grid: [" << dim_grid(0)
            << ", " << dim_grid(1) << ", " << dim_grid(2)
            << "], total vol: " << (vol_flows + vol_parts + vol_mesh_flows)
            << ", " << vol_parts << ", " << vol_flows << ", " << vol_mesh_flows
            << "]" << std::endl;
}

FluidScriptingController::FluidTiming operator+(
    const FluidScriptingController::FluidTiming& lhs,
    const FluidScriptingController::FluidTiming& rhs) {
  FluidScriptingController::FluidTiming sum;
  sum.processparts = lhs.processparts + rhs.processparts;
  sum.applypressure = lhs.applypressure + rhs.applypressure;
  sum.stress = lhs.stress + rhs.stress;
  sum.datastructure = lhs.datastructure + rhs.datastructure;
  sum.solidproj = lhs.solidproj + rhs.solidproj;
  sum.plasticity = lhs.plasticity + rhs.plasticity;
  sum.weight = lhs.weight + rhs.weight;
  sum.advect = lhs.advect + rhs.advect;
  sum.mapp2g = lhs.mapp2g + rhs.mapp2g;
  sum.addforce = lhs.addforce + rhs.addforce;
  sum.computephi = lhs.computephi + rhs.computephi;
  sum.solvedrag = lhs.solvedrag + rhs.solvedrag;
  sum.solvevisc = lhs.solvevisc + rhs.solvevisc;
  sum.solvepressure = lhs.solvepressure + rhs.solvepressure;
  sum.constraintvel = lhs.constraintvel + rhs.constraintvel;
  sum.mapg2p = lhs.mapg2p + rhs.mapg2p;
  sum.syncflow = lhs.syncflow + rhs.syncflow;

  return sum;
}

template <typename StreamT>
void FluidScriptingController::print(
    const FluidScriptingController::FluidTiming& timings) const {
  StreamT(g_log, "Fluid Timing")
      << " PP: " << timings.processparts << " DS: " << timings.datastructure
      << " CW: " << timings.weight << " CP: " << timings.computephi
      << " P2G: " << timings.mapp2g << " AF: " << timings.addforce
      << " SV: " << timings.solvevisc << " SP: " << timings.solvepressure
      << " AP: " << timings.applypressure << " SF: " << timings.syncflow
      << " CV: " << timings.constraintvel << " G2P: " << timings.mapg2p
      << " SLP: " << timings.solidproj << " ADV: " << timings.advect
      << " UP: " << timings.plasticity;

  StreamT(g_log, "Fluid Timing") << "Total: " << timings.sum();
}

FluidScriptingController::FluidTiming operator/(
    const FluidScriptingController::FluidTiming& lhs, const Scalar rhs) {
  FluidScriptingController::FluidTiming avg = lhs;
  avg.processparts /= rhs;
  avg.applypressure /= rhs;
  avg.stress /= rhs;
  avg.datastructure /= rhs;
  avg.solidproj /= rhs;
  avg.plasticity /= rhs;
  avg.weight /= rhs;
  avg.advect /= rhs;
  avg.mapp2g /= rhs;
  avg.addforce /= rhs;
  avg.computephi /= rhs;
  avg.solvedrag /= rhs;
  avg.solvevisc /= rhs;
  avg.solvepressure /= rhs;
  avg.constraintvel /= rhs;
  avg.mapg2p /= rhs;
  avg.syncflow /= rhs;

  return avg;
}

}  // namespace strandsim
