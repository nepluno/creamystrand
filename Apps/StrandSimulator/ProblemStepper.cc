/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "ProblemStepper.hh"

#include <time.h>

#include "../../StrandSim/Dynamic/ImplicitStepper.hh"
#include "../../StrandSim/Dynamic/TriangularMeshFlow.hh"

void dump_data_subprog(DumpData* data);

ProblemStepper::ProblemStepper(const std::string& name, const std::string& desc)
    : m_problemName(name),
      m_problemDesc(desc),
      m_t(0.)  // start time
      ,
      m_isSimulated(true),
      m_stepper(NULL),
      m_num_components(1),
      m_cfl_number(1.0) {}

ProblemStepper::~ProblemStepper() {
  if (m_stepper != NULL) {
    delete m_stepper;
    m_stepper = NULL;
  }

  for (std::vector<strandsim::ElasticStrand*>::size_type i = 0;
       i < m_strands.size(); ++i) {
    if (m_strands[i] != NULL) {
      delete m_strands[i];
      m_strands[i] = NULL;
    }
  }

  for (std::vector<RodData*>::size_type i = 0; i < m_rodDatum.size(); ++i) {
    delete m_rodDatum[i];
    m_rodDatum[i] = NULL;
  }

  for (std::vector<strandsim::TriangleMeshRenderer*>::size_type i = 0;
       i < m_mesh_renderers.size(); ++i) {
    if (m_mesh_renderers[i] != NULL) {
      delete m_mesh_renderers[i];
      m_mesh_renderers[i] = NULL;
    }
  }

  for (std::vector<strandsim::FluidsRenderer*>::size_type i = 0;
       i < m_fluids_renderers.size(); ++i) {
    if (m_fluids_renderers[i] != NULL) {
      delete m_fluids_renderers[i];
      m_fluids_renderers[i] = NULL;
    }
  }

  // Delete these too?
  // std::vector<ConstraintScriptingController*>
  // m_constraintScripting_controllers;

  for (std::vector<DOFScriptingController*>::size_type i = 0;
       i < m_dof_scripting_controllers.size(); ++i) {
    if (m_dof_scripting_controllers[i] != NULL) {
      delete m_dof_scripting_controllers[i];
      m_dof_scripting_controllers[i] = NULL;
    }
  }
}

void ProblemStepper::setup(int& current_frame, int& current_check_point) {
  // std::cout << "Set up!\n";
  setSimulationParameters();
  setupStrands();
  setupMeshes();

  if (m_isSimulated)
    m_stepper = new StrandImplicitManager(
        m_strands, m_collision_free, m_meshScripting_controllers,
        m_fluidScripting_controllers, m_constraintScripting_controllers, m_t,
        m_dt, m_simulation_params, this);

  setupAfterInit(current_frame, current_check_point);
  //    for( auto rod_id = m_strands.begin(); rod_id != m_strands.end();
  //    ++rod_id)
  //    {
  //        std::cout << "Initial vertices for rod " <<
  //        (*rod_id)->getGlobalIndex() << " : " << *(*rod_id) << '\n';
  //    }
}
namespace strandsim {}
bool ProblemStepper::step() {
  // by default use an empirical CFL number 1.0 for shear-dependent fluid
  const Scalar cfl = m_stepper->getCFL();
  const Scalar dt_cfl = cfl * m_cfl_number;

  int num_substeps = std::max(1, (int)ceil(m_dt / dt_cfl));

  const Scalar dt_total = m_dt / (Scalar)num_substeps;

  for (int i_step = 0; i_step < num_substeps; ++i_step) {
    std::cout << "\n# Stepping [" << i_step << "/" << num_substeps
              << "] @ time: " << m_t << " with " << dt_total
              << " (max: " << dt_cfl << ")" << std::endl;
    m_t += dt_total;

    if (!executeScript(i_step, dt_total)) return false;

    if (m_isSimulated) m_stepper->execute(num_substeps, i_step, dt_total);
  }

  return true;
}

void ProblemStepper::executeCallback() {
  for (auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr) {
    (*rd_itr)->update();
  }
}

void ProblemStepper::render() {
  for (auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr) {
    (*rd_itr)->render();
  }

  for (auto m_itr = m_mesh_renderers.begin(); m_itr != m_mesh_renderers.end();
       ++m_itr) {
    (*m_itr)->render();
  }

  for (auto m_itr = m_fluids_renderers.begin();
       m_itr != m_fluids_renderers.end(); ++m_itr) {
    (*m_itr)->render();
  }
}

void ProblemStepper::getCenter(Vec3d& center) {
  int render_count = 0;
  for (auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr) {
    center += (*rd_itr)->calculateObjectCenter();
    ++render_count;
  }

  for (auto m_itr = m_mesh_renderers.begin(); m_itr != m_mesh_renderers.end();
       ++m_itr) {
    center += (*m_itr)->calculateObjectCenter();
    ++render_count;
  }

  if (render_count > 0) center /= render_count;
}

void ProblemStepper::getRadius(Scalar& radius, const Vec3d& simCenter) {
  for (auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr) {
    const Vec3x center = (*rd_itr)->calculateObjectCenter();
    const Scalar r = (*rd_itr)->calculateObjectBoundingRadius(center);
    radius = std::max(radius, r + (center - simCenter).norm());
  }

  for (auto m_itr = m_mesh_renderers.begin(); m_itr != m_mesh_renderers.end();
       ++m_itr) {
    const Vec3x center = (*m_itr)->calculateObjectCenter();
    const Scalar r = (*m_itr)->calculateObjectBoundingRadius(center);
    radius = std::max(radius, r + (center - simCenter).norm());
  }
}

void ProblemStepper::findOrthogonal(Vec3x& v, const Vec3x& u) {
  assert(u.norm() != 0);

  v.setZero();
  int max = 0;
  for (int i = 0; i < u.size(); ++i) {
    if (u[i] == 0) {
      v[i] = 1;
      return;
    }
    if (fabs(u[i]) > fabs(u[max])) max = i;
  }

  int idx = (max + 1) % u.size();
  v[idx] = u[max];
  v[max] = -u[idx];
  v.normalize();

  assert(std::abs(u.dot(v)) < 1e-10);
  // assert(u.dot(v)== 0.0);
}

void ProblemStepper::transformTriangleObject(TriangularMesh& triMesh,
                                             const Mat3x& transformation,
                                             const Vec3x& center,
                                             const Vec3x& translate) {
  for (unsigned i = 0; i < triMesh.nv(); ++i) {
    Vec3x vert = triMesh.getVertex(i);
    Vec3x vertNext = transformation * (vert - center) + center + translate;
    triMesh.setVertex(i, vertNext);
    triMesh.setDisplacement(i, (vertNext - vert));
  }
}

void ProblemStepper::freezeTriangleObject(TriangularMesh& triMesh) {
  Vec3x disp(0., 0., 0.);
  for (unsigned i = 0; i < triMesh.nv(); ++i) {
    triMesh.setDisplacement(i, disp);
  }
}

void ProblemStepper::transformRodRootVtx(RodData& rd,
                                         const Mat3x& transformation,
                                         const Vec3x& center,
                                         const Vec3x& translate,
                                         const Vec3x& scale, int vtx) {
  if (m_simulation_params.m_useSoftAttachConstraints) {
    Vec3x vert = rd.getDofController().getVertexGoal(
        vtx);  // rd.getStrand().getVertex( vtx );
    Vec3x vertNext =
        transformation * Vec3x((vert - center).array() * scale.array()) +
        center + translate;
    rd.getDofController().setVertexGoal(vtx, vertNext);
  } else {
    Vec3x vert = rd.getStrand().getVertex(vtx);
    Vec3x vertNext =
        transformation * Vec3x((vert - center).array() * scale.array()) +
        center + translate;
    rd.getDofController().setVertexDisplacement(
        vtx, (vertNext - vert) / (Scalar)m_simulation_params.m_subSteps);
  }
}

void ProblemStepper::transformTriangleObject(
    TriangularMesh& triMesh, const Eigen::Quaternion<double>& transformation,
    const Vec3x& center, const Vec3x& translate) {
  for (unsigned i = 0; i < triMesh.nv(); ++i) {
    Vec3x vert = triMesh.getVertex(i);
    Vec3x vertNext = transformation * (vert - center) + center + translate;
    triMesh.setVertex(i, vertNext);
    triMesh.setDisplacement(i, (vertNext - vert));
  }
}

void ProblemStepper::transformRodRootVtx(
    RodData& rd, const Eigen::Quaternion<double>& transformation,
    const Vec3x& center, const Vec3x& translate, const Vec3x& scale, int vtx) {
  if (m_simulation_params.m_useSoftAttachConstraints) {
    Vec3x vert = rd.getDofController().getVertexGoal(vtx);
    Vec3x vertNext =
        transformation * Vec3x((vert - center).array() * scale.array()) +
        center + translate;
    rd.getDofController().setVertexGoal(vtx, vertNext);
  } else {
    Vec3x vert = rd.getStrand().getVertex(vtx);
    Vec3x vertNext =
        transformation * Vec3x((vert - center).array() * scale.array()) +
        center + translate;
    rd.getDofController().setVertexDisplacement(
        vtx, (vertNext - vert) / (Scalar)m_simulation_params.m_subSteps);
  }
}

void ProblemStepper::setStatLogging(int statLogging) {
  m_simulation_params.m_statGathering = statLogging;
}

void ProblemStepper::PrintAdditionalSettings(const std::string& szfn_dir) {}

void ProblemStepper::dumpFluid(DumpData* data, std::string outputdirectory,
                               int current_frame, int file_width) const {
  int mesh_num = 0;
  for (auto m_itr = m_fluidScripting_controllers.begin();
       m_itr != m_fluidScripting_controllers.end(); ++m_itr, ++mesh_num) {
    // new obj per mesh
    std::stringstream name;
    name << std::setfill('0');
    name << outputdirectory << "/fluid" << mesh_num << "_"
         << std::setw(file_width) << current_frame << ".ply";

    data->fluid_fn = name.str().c_str();

    auto controller = *m_itr;
    const int numv = controller->numParticles();

    data->fluid_particles.resize(numv);
    data->fluid_radius.resize(numv);
    data->fluid_colors.resize(numv);

    for (size_t v = 0; v < numv; ++v) {
      data->fluid_particles[v] = controller->getParticles(v);
      data->fluid_radius[v] = controller->getRadius(v);
      data->fluid_colors[v] = controller->getComponents(v);
    }
  }
}

void ProblemStepper::dumpMesh(DumpData* data, std::string outputdirectory,
                              int current_frame, int file_width) const {
  int mesh_num = 0;
  data->mesh_indices.resize(m_meshScripting_controllers.size());
  data->mesh_vertices.resize(m_meshScripting_controllers.size());
  data->mesh_fn.resize(m_meshScripting_controllers.size());
  data->mesh_flow_height.resize(m_meshScripting_controllers.size());
  data->mesh_colors.resize(m_meshScripting_controllers.size());

  for (auto m_itr = m_meshScripting_controllers.begin();
       m_itr != m_meshScripting_controllers.end(); ++m_itr, ++mesh_num) {
    // new obj per mesh
    std::stringstream name;
    name << std::setfill('0');
    name << outputdirectory << "/mesh" << mesh_num << "_"
         << std::setw(file_width) << current_frame << ".ply";
    data->mesh_fn[mesh_num] = name.str().c_str();

    auto mesh = (*m_itr)->getCurrentMesh();
    auto flow = (*m_itr)->getMeshFlow();

    const int nv = mesh->nv();
    const int nf = mesh->nf();
    data->mesh_vertices[mesh_num].resize(mesh->nv());
    data->mesh_flow_height[mesh_num].resize(mesh->nv());
    data->mesh_indices[mesh_num].resize(mesh->nf());
    data->mesh_colors[mesh_num].resize(mesh->nv());

    for (int v = 0; v < nv; ++v) {
      data->mesh_vertices[mesh_num][v] = mesh->getVertex(v);
    }

    if (flow) {
      const VecXx& flow_height = flow->getCurrentFlowHeight();
      const VecXx& flow_colors = flow->getCurrentComponents();

      for (int v = 0; v < nv; ++v) {
        data->mesh_flow_height[mesh_num][v] = flow_height[v];
        data->mesh_colors[mesh_num][v] =
            flow_colors.segment(v * m_num_components, m_num_components);
      }
    } else {
      data->mesh_flow_height[mesh_num].assign(nv, 0.0);
      data->mesh_colors[mesh_num].assign(nv, VecXx::Zero(m_num_components));
    }

    for (int f = 0; f < nf; ++f) {
      data->mesh_indices[mesh_num][f] =
          Vec3i(mesh->getFace(f).idx[0], mesh->getFace(f).idx[1],
                mesh->getFace(f).idx[2]);
    }
  }
}

void ProblemStepper::setOutputDirectory(const std::string& dir) {
  if (!m_stepper) return;
  m_stepper->setOutputDirectory(dir);
}

void ProblemStepper::dumpBinaryCheckpoint(std::string outputdirectory,
                                          int current_frame,
                                          int current_checkpoint,
                                          int file_width) const {}

void ProblemStepper::dumpRods(DumpData* data, std::string outputdirectory,
                              int current_frame, int file_width) const {
  std::stringstream name;
  name << std::setfill('0');
  name << outputdirectory << "/rods_" << std::setw(file_width) << current_frame
       << ".ply";

  data->rod_fn = name.str().c_str();

  // we add extra 2x2 vertices before first and after last, to have a nice cap
  // for the flow during reconstruction

  int num_verts = 0;
  for (auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr) {
    num_verts += (*rd_itr)->getStrand().getNumVertices() + 4;
  }

  data->rod_radius.reserve(num_verts);
  data->rod_indices.reserve(num_verts);
  data->rod_vertices.reserve(num_verts);
  data->rod_flow_height.reserve(num_verts);
  data->rod_groups.reserve(num_verts);
  data->rod_colors.reserve(num_verts);
  data->rod_start_vert_indices.reserve(m_rodDatum.size());

  int rod_num = 0;
  for (auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end();
       ++rd_itr, ++rod_num) {
    const ElasticStrand& strand = (*rd_itr)->getStrand();
    // ignore bad strands
    if (strand.getNumVertices() < 2) continue;

    data->rod_start_vert_indices.push_back((int)(data->rod_vertices.size()));

    Scalar h0 = strand.getCurrentFlowHeight(0);
    Scalar r0 = sqrt(strand.getRadiusA(0) * strand.getRadiusB(0));
    Scalar hr0 = h0 + r0;
    Scalar r_rs0 = pow(strand.getReservoir()(0) * 0.75 / M_PI, 1.0 / 3.0);
    Scalar e0 = sqrt(std::max(r_rs0 * r_rs0 - hr0 * hr0, hr0 * hr0));
    Vec3x t0 = -strand.getEdgeVector(0).normalized();
    Vec3x v00 = strand.getVertex(0) + e0 * t0;
    Vec3x v01 = v00 + r_rs0 * t0;
    VecXx c0 =
        strand.getStepper()->flowComponents().segment(0, m_num_components);

    data->rod_vertices.push_back(Vec4x(v01(0), v01(1), v01(2), 0.0));
    data->rod_indices.push_back(rod_num);
    data->rod_radius.push_back(Vec2x::Zero());
    data->rod_flow_height.push_back(0.0);
    data->rod_groups.push_back(strand.getParameters().getParamIndex());
    data->rod_actual.push_back(0);
    data->rod_colors.push_back(c0);

    data->rod_vertices.push_back(Vec4x(v00(0), v00(1), v00(2), 0.0));
    data->rod_indices.push_back(rod_num);
    data->rod_radius.push_back(Vec2x::Zero());
    data->rod_flow_height.push_back(r_rs0);
    data->rod_groups.push_back(strand.getParameters().getParamIndex());
    data->rod_actual.push_back(0);
    data->rod_colors.push_back(c0);

    const int n = strand.getNumVertices();
    for (int j = 0; j < n; ++j) {
      Vec3x vv = strand.getVertex(j);
      Scalar vt = (j == n - 1) ? 0.0 : strand.getTheta(j);

      data->rod_vertices.push_back(Vec4x(vv(0), vv(1), vv(2), vt));
      data->rod_indices.push_back(rod_num);
      data->rod_radius.push_back(
          Vec2x(strand.getRadiusA(j), strand.getRadiusB(j)));
      data->rod_flow_height.push_back(strand.getCurrentFlowHeight(j));
      data->rod_groups.push_back(strand.getParameters().getParamIndex());
      data->rod_actual.push_back(1);
      data->rod_colors.push_back(
          VecXx(strand.getStepper()->flowComponents().segment(
              j * m_num_components, m_num_components)));
    }

    Scalar h1 = strand.getCurrentFlowHeight(n - 1);
    Scalar r1 = sqrt(strand.getRadiusA(n - 1) * strand.getRadiusB(n - 1));
    Scalar hr1 = h1 + r1;
    Scalar r_rs1 = pow(strand.getReservoir()(1) * 0.75 / M_PI, 1.0 / 3.0);
    Scalar e1 = sqrt(std::max(r_rs1 * r_rs1 - hr1 * hr1, hr1 * hr1));
    Vec3x t1 = strand.getEdgeVector(n - 2).normalized();
    Vec3x v10 = strand.getVertex(n - 1) + e1 * t1;
    Vec3x v11 = v10 + r_rs1 * t1;
    VecXx c1 = strand.getStepper()->flowComponents().segment(
        (n - 1) * m_num_components, m_num_components);

    data->rod_vertices.push_back(Vec4x(v10(0), v10(1), v10(2), 0.0));
    data->rod_indices.push_back(rod_num);
    data->rod_radius.push_back(Vec2x::Zero());
    data->rod_flow_height.push_back(r_rs1);
    data->rod_groups.push_back(strand.getParameters().getParamIndex());
    data->rod_actual.push_back(0);
    data->rod_colors.push_back(c1);

    data->rod_vertices.push_back(Vec4x(v11(0), v11(1), v11(2), 0.0));
    data->rod_indices.push_back(rod_num);
    data->rod_radius.push_back(Vec2x::Zero());
    data->rod_flow_height.push_back(0.0);
    data->rod_groups.push_back(strand.getParameters().getParamIndex());
    data->rod_actual.push_back(0);
    data->rod_colors.push_back(c1);
  }
}

void ProblemStepper::dumpData(std::string outputdirectory, int current_frame,
                              int file_width) const {
  DumpData* data = new DumpData;
  dumpFluid(data, outputdirectory, current_frame, file_width);
  dumpMesh(data, outputdirectory, current_frame, file_width);
  dumpRods(data, outputdirectory, current_frame, file_width);

  data->current_frame = current_frame;
  data->num_components = m_num_components;

  std::thread t(std::bind(dump_data_subprog, data));
  t.detach();
}

void ProblemStepper::projectConstraint(
    const std::vector<ImplicitStepper*>& steppers) {}

void dump_data_subprog(DumpData* data) {
  std::ofstream os_fluid(data->fluid_fn.c_str());

  // header
  os_fluid << "ply" << std::endl
           << "format ascii 1.0" << std::endl
           << "comment created by ADONIS" << std::endl;
  int num_parts = (int)data->fluid_particles.size();
  os_fluid << "element vertex " << num_parts << std::endl;
  os_fluid << "property float x" << std::endl
           << "property float y" << std::endl
           << "property float z" << std::endl
           << "property float pscale" << std::endl;

  for (int i = 0; i < data->num_components; ++i) {
    os_fluid << "property float c" << i << std::endl;
  }

  os_fluid << "element face 0" << std::endl
           << "property list int int vertex_indices" << std::endl
           << "end_header " << std::endl;

  const int numv = (int)data->fluid_particles.size();
  for (int v = 0; v < numv; ++v) {
    const Vec3x& p = data->fluid_particles[v];
    os_fluid << p(0) << " " << p(1) << " " << p(2) << " "
             << data->fluid_radius[v] << " ";

    for (int j = 0; j < data->num_components; ++j) {
      os_fluid << data->fluid_colors[v][j];
      if (j != data->num_components - 1) {
        os_fluid << " ";
      }
    }

    os_fluid << std::endl;
  }

  os_fluid.flush();
  os_fluid.close();

  const int num_meshes = (int)data->mesh_fn.size();
  for (int mesh_num = 0; mesh_num < num_meshes; ++mesh_num) {
    std::ofstream os(data->mesh_fn[mesh_num].c_str());

    os << "ply" << std::endl
       << "format ascii 1.0" << std::endl
       << "comment created by ADONIS" << std::endl;
    const int num_verts = (int)data->mesh_vertices[mesh_num].size();
    const int num_faces = (int)data->mesh_indices[mesh_num].size();

    os << "element vertex " << num_verts << std::endl;
    os << "property float x" << std::endl
       << "property float y" << std::endl
       << "property float z" << std::endl
       << "property float h" << std::endl;

    for (int j = 0; j < data->num_components; ++j) {
      os << "property float c" << j << std::endl;
    }

    os << "element face " << num_faces << std::endl
       << "property list int int vertex_indices" << std::endl
       << "end_header " << std::endl;

    for (int v = 0; v < num_verts; ++v) {
      os << data->mesh_vertices[mesh_num][v](0) << " "
         << data->mesh_vertices[mesh_num][v](1) << " "
         << data->mesh_vertices[mesh_num][v](2) << " "
         << data->mesh_flow_height[mesh_num][v] << " ";

      for (int j = 0; j < data->num_components; ++j) {
        os << data->mesh_colors[mesh_num][v][j];

        if (j != data->num_components - 1) {
          os << " ";
        }
      }

      os << std::endl;
    }

    for (int f = 0; f < num_faces; ++f) {
      os << "3 " << data->mesh_indices[mesh_num][f](0) << " "
         << data->mesh_indices[mesh_num][f](1) << " "
         << data->mesh_indices[mesh_num][f](2) << std::endl;
    }

    os.flush();
    os.close();
  }

  std::ofstream os_rod(data->rod_fn.c_str());

  // header
  os_rod << "ply" << std::endl
         << "format ascii 1.0" << std::endl
         << "comment created by ADONIS" << std::endl;
  int num_verts = (int)data->rod_vertices.size();
  os_rod << "element vertex " << num_verts << std::endl;
  os_rod << "property float x" << std::endl
         << "property float y" << std::endl
         << "property float z" << std::endl
         << "property float theta" << std::endl
         << "property int segment" << std::endl
         << "property float ra" << std::endl
         << "property float rb" << std::endl
         << "property float ha" << std::endl
         << "property float hb" << std::endl
         << "property int group" << std::endl
         << "property int actual" << std::endl;

  for (int j = 0; j < data->num_components; ++j) {
    os_rod << "property float c" << j << std::endl;
  }

  const int num_rods = (int)data->rod_start_vert_indices.size();

  os_rod << "element face " << num_rods << std::endl;
  os_rod << "property list int int vertex_indices" << std::endl
         << "end_header " << std::endl;

  for (int i = 0; i < num_verts; ++i) {
    os_rod << data->rod_vertices[i](0) << " " << data->rod_vertices[i](1) << " "
           << data->rod_vertices[i](2) << " " << data->rod_vertices[i](3) << " "
           << data->rod_indices[i] << " " << data->rod_radius[i](0) << " "
           << data->rod_radius[i](1) << " "
           << (data->rod_radius[i](0) + data->rod_flow_height[i]) << " "
           << (data->rod_radius[i](1) + data->rod_flow_height[i]) << " "
           << data->rod_groups[i] << " " << data->rod_actual[i] << " ";

    for (int j = 0; j < data->num_components; ++j) {
      os_rod << data->rod_colors[i][j];

      if (j != data->num_components - 1) {
        os_rod << " ";
      }
    }

    os_rod << std::endl;
  }

  for (int i = 0; i < num_rods; ++i) {
    const int num_vert_in_rod =
        (i == num_rods - 1)
            ? ((int)data->rod_vertices.size() - data->rod_start_vert_indices[i])
            : (data->rod_start_vert_indices[i + 1] -
               data->rod_start_vert_indices[i]);
    os_rod << num_vert_in_rod << " ";

    for (int j = 0; j < num_vert_in_rod; ++j) {
      os_rod << (data->rod_start_vert_indices[i] + j) << " ";
    }

    os_rod << "0" << std::endl;
  }

  os_rod.flush();
  os_rod.close();

  std::cout << "[Frame " << data->current_frame << " written]" << std::endl;
  delete data;
}
