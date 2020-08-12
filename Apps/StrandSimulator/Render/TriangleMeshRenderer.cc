/**
 * \copyright 2010 Breannan Smith
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "TriangleMeshRenderer.hh"

#include "../../StrandSim/Dynamic/FluidScriptingController.hh"
#include "../../StrandSim/Dynamic/TriangularMeshFlow.hh"
#include "../../StrandSim/Render/Color.hh"
#include "../../StrandSim/Render/OpenGLDecl.hh"
#include "../../StrandSim/Utils/ThreadUtils.hh"

namespace strandsim {

TriangleMeshRenderer::TriangleMeshRenderer(TriangularMesh& mesh,
                                           int num_components)
    : m_mesh(mesh), m_mode(FLAT), m_controller(NULL) {
  auto idx2rgb = [&](int i) -> Vec3x {
    Scalar hh = (Scalar)i / (Scalar)std::max(1, num_components - 1) * 5.0;
    if (hh == 6.0) hh = 0.0;

    int fi = (int)floor(hh);
    Scalar ff = hh - floor(hh);

    switch (fi) {
      case 0:
        return Vec3x(1.0, ff, 0.0);
      case 1:
        return Vec3x(1.0 - ff, 1.0, 0.0);
      case 2:
        return Vec3x(0.0, 1.0, ff);
      case 3:
        return Vec3x(0.0, 1.0 - ff, 1.0);
      case 4:
        return Vec3x(ff, 0.0, 1.0);
      case 5:
      default:
        return Vec3x(1.0, 0.0, 1.0 - ff);
    }
  };

  m_component_colors.resize(num_components);
  for (int i = 0; i < num_components; ++i) {
    m_component_colors[i] = idx2rgb(i);
  }
}

void TriangleMeshRenderer::setFluidController(
    const std::shared_ptr<FluidScriptingController>& controller) {
  m_controller = controller;
}

void TriangleMeshRenderer::render() {
  LockGuard lock(m_mesh.getGeometryMutex());

  glEnable(GL_BLEND);
  glDisable(GL_DEPTH_TEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  if (m_mode == FLAT) {
    glEnable(GL_LIGHTING);

    GLfloat gray[] = {0.8, 0.8, 0.8, 0.25};
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, gray);

    // Render all faces
    glBegin(GL_TRIANGLES);
    // OpenGL::color(Color(255,0,0));

    //    for( TriangleMesh::face_iter fit = m_mesh.faces_begin(); fit !=
    //    m_mesh.faces_end(); ++fit )
    //    {
    //      std::vector<Vec3d> v;
    //      for( TriangleMesh::FaceVertexIter fvit = m_mesh.fv_iter(*fit); fvit;
    //      ++fvit )
    //      {
    //        v.push_back(m_mesh.getVertex(*fvit));
    //      }

    int fidx = 0;
    for (auto fit = m_mesh.getFaces().begin(); fit != m_mesh.getFaces().end();
         ++fit) {
      std::vector<Vec3d> v;
      for (int i = 0; i < 3; ++i) v.push_back(m_mesh.getVertex((*fit).idx[i]));

      // Compute a normal for the face
      Vec3d n = (v[1] - v[0]).cross(v[2] - v[0]).normalized();

      glNormal3f(n.x(), n.y(), n.z());
      glVertex3f(v[0].x(), v[0].y(), v[0].z());
      glVertex3f(v[1].x(), v[1].y(), v[1].z());
      glVertex3f(v[2].x(), v[2].y(), v[2].z());

      ++fidx;
    }
    glEnd();

    if (m_mesh.m_associatedFlow) {
      //                glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
      const VecXx& flow_height =
          m_mesh.m_associatedFlow->getCurrentFlowHeight();
      const VecXx& flow_vel = m_mesh.m_associatedFlow->getCurrentFlowVelocity();

      // Render all faces
      glBegin(GL_TRIANGLES);

      int fidx = 0;
      for (auto fit = m_mesh.getFaces().begin(); fit != m_mesh.getFaces().end();
           ++fit) {
        std::vector<Vec3d> v;

        //                    double accum = 0.0;

        for (int i = 0; i < 3; ++i) {
          int vidx = (*fit).idx[i];
          Vec3d nv = m_mesh.m_vertex_normals[vidx];
          v.push_back(m_mesh.getVertex(vidx) + nv * flow_height(vidx));
          //                        accum += flow_height(vidx);
        }

        // Compute a normal for the face
        for (int r = 0; r < 3; ++r) {
          const int idx_vert = (*fit).idx[r];
          Vec3d n = m_mesh.m_vertex_normals[idx_vert];
          glNormal3dv(n.data());

          if (!m_controller) {
            GLfloat blue[] = {0.0f, 0.0f, 0.8f,
                              (float)flow_height(idx_vert) / 0.1f};
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, blue);
          } else {
            const VecXx& c = m_controller->getComponents(idx_vert);
            Vec3x color = Vec3x::Zero();
            for (int j = 0; j < c.size(); ++j) {
              color += m_component_colors[j] * c(j);
            }
            GLfloat blue[] = {(GLfloat)color(0), (GLfloat)color(1),
                              (GLfloat)color(2),
                              (GLfloat)flow_height(idx_vert) / 0.1f};
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, blue);
          }
          glVertex3dv(v[r].data());
        }

        ++fidx;
      }
      glEnd();
      //                glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

      GLfloat blue[] = {0.0, 0.0, 0.8, .25f};
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, blue);
      glBegin(GL_LINES);
      const int num_faces = m_mesh.nf();
      for (int i = 0; i < num_faces; ++i) {
        const TriangularFace& f = m_mesh.getFace(i);
        const Vec3x center =
            (m_mesh.getVertex(f(0)) +
             m_mesh.m_vertex_normals[f(0)] * flow_height(f(0)) +
             m_mesh.getVertex(f(1)) +
             m_mesh.m_vertex_normals[f(1)] * flow_height(f(1)) +
             m_mesh.getVertex(f(2)) +
             m_mesh.m_vertex_normals[f(2)] * flow_height(f(2))) /
            3.0;
        const Vec3x vel =
            m_mesh.m_associatedFlow->getCurrentFlowVelocity().segment<3>(i * 3);
        const Vec3x target = center + vel * 0.01;

        glVertex3dv(center.data());
        glVertex3dv(target.data());
      }

      glEnd();
    }

    glDisable(GL_LIGHTING);
  } else if (m_mode == DBG) {
    glDisable(GL_LIGHTING);

    //    // Render all edges
    //    glLineWidth(2);
    //    glBegin(GL_LINES);
    //    OpenGL::color(Color(0,0,0));
    //    for( TriangleMesh::edge_iter eit = m_mesh.edges_begin(); eit !=
    //    m_mesh.edges_end(); ++eit )
    //    {
    //      OpenGL::vertex(m_mesh.getVertex(m_mesh.fromVertex(*eit)));
    //      OpenGL::vertex(m_mesh.getVertex(m_mesh.toVertex(*eit)));
    //    }
    //    glEnd();

    // Render all faces
    glBegin(GL_TRIANGLES);
    OpenGL::color(Color(255, 0, 0));
    //    for( TriangleMesh::face_iter fit = m_mesh.faces_begin(); fit !=
    //    m_mesh.faces_end(); ++fit )
    //    {
    //      for( TriangleMesh::FaceVertexIter fvit = m_mesh.fv_iter(*fit); fvit;
    //      ++fvit )
    //      {
    for (auto fit = m_mesh.getFaces().begin(); fit != m_mesh.getFaces().end();
         ++fit) {
      for (int i = 0; i < 3; ++i) {
        OpenGL::vertex(m_mesh.getVertex((*fit).idx[i]));
        //        OpenGL::vertex(m_mesh.getVertex(*fvit));
      }
    }
    glEnd();

    // Render all vertices
    glPointSize(5);
    glBegin(GL_POINTS);
    OpenGL::color(Color(0, 0, 0));
    //    for( TriangleMesh::vertex_iter vit = m_mesh.vertices_begin(); vit !=
    //    m_mesh.vertices_end(); ++vit ) OpenGL::vertex(m_mesh.getVertex(*vit));
    for (unsigned i = 0; i < m_mesh.nv(); ++i) {
      OpenGL::vertex(m_mesh.getVertex(i));
    }

    glEnd();

    glEnable(GL_LIGHTING);
  }
  glDisable(GL_CULL_FACE);
  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
}

Vec3d TriangleMeshRenderer::calculateObjectCenter() {
  Vec3d center(0.0, 0.0, 0.0);

  //  for( TriangleMesh::vertex_iter vit = m_mesh.vertices_begin(); vit !=
  //  m_mesh.vertices_end(); ++vit )
  //  {
  //    center += m_mesh.getVertex(*vit);
  //  }
  for (unsigned i = 0; i < m_mesh.nv(); ++i) {
    center += m_mesh.getVertex(i);
  }

  if (m_mesh.nv() != 0) center /= ((double)m_mesh.nv());

  return center;
}

double TriangleMeshRenderer::calculateObjectBoundingRadius(
    const Vec3d& center) {
  Scalar radius = 0.0;

  //  for( TriangleMesh::vertex_iter vit = m_mesh.vertices_begin(); vit !=
  //  m_mesh.vertices_end(); ++vit )
  //  {
  //    radius = std::max(radius, (m_mesh.getVertex(*vit) - center).norm());
  //  }

  for (unsigned i = 0; i < m_mesh.nv(); ++i) {
    radius = std::max(radius, (m_mesh.getVertex(i) - center).norm());
  }

  return radius;
}

}  // namespace strandsim
