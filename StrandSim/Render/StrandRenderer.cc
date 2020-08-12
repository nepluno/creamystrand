/**
 * \copyright 2011 Susan Howard, Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "StrandRenderer.hh"

#include "../Core/ElasticStrand.hh"
#include "../Core/ElasticStrandUtils.hh"
#include "../Core/StrandState.hh"
#include "../Dynamic/FluidScriptingController.hh"
#include "../Dynamic/ImplicitStepper.hh"
#include "../Dynamic/StrandImplicitManager.hh"
#include "../Forces/BendingForce.hh"
#include "../Forces/ForceAccumulator.hh"
#include "../Forces/GravitationForce.hh"
#include "../Forces/StretchingForce.hh"
#include "../Forces/TwistingForce.hh"
#include "../Utils/MathUtilities.hh"

namespace strandsim {

std::set<int> StrandRenderer::s_debugDisplaySet;

StrandRenderer::StrandRenderer(const ElasticStrand& strand,
                               const std::vector<float>& currentVertices,
                               const std::vector<float>& currentMaterialFrames,
                               MutexType& geometryMutex, int num_components)
    : m_strand(strand),                                //
      m_currentVertices(currentVertices),              //
      m_currentMaterialFrames(currentMaterialFrames),  //
      m_selected(false),                               //
      m_simulated(true),                               //
      m_immune(false),
      m_drawStrand(true),
      m_drawFlow(true),           //
      m_drawRootMaterial(false),  //
      m_drawAnchors(false),       //
      m_drawMode(LINES),  // QUADS_SHADED ), //QUADS_SHADED), //LINES ), //
      m_drawForce(NONE),  //
      m_drawnVolume(NO_VOLUME),  //
      m_geometryChanged(true),   //
      m_flowChanged(true),
      m_colorChanged(true),            //
      m_geometryMutex(geometryMutex),  //
      m_drawDebug(false),
      m_lineWidth(1) {
  m_transformationMatrix.setIdentity();

  m_palette["first material"] = Color(255, 125, 75);
  m_palette["second material"] = Color(55, 125, 255);
  m_palette["first reference frame"] = Color(200, 200, 200);
  m_palette["second reference frame"] = Color(50, 50, 50);
  m_palette["velocity"] = Color(0, 0, 255);
  m_palette["force"] = Color(255, 0, 0);
  m_palette["simple root"] = Color(20, 25, 13);  // Color( 0, 0, 0 );
  m_palette["simple tip"] =
      Color(169, 99, 49);  // Color(139, 69, 19); //Color( 155, 200, 100 );
  m_palette["default simple tip"] = m_palette["simple tip"];
  m_palette["external collisions volume"] = Color(155, 200, 100, 128);
  m_palette["self collisions volume"] = Color(125, 195, 255);
  m_palette["water"] = Color(64, 164, 223);
  m_palette["random material"] =
      Color(((float)(rand() % 255)) / 255., ((float)(rand() % 255)) / 255.,
            ((float)(rand() % 255)) / 255.);

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

void StrandRenderer::render() {
  glPushAttrib(GL_COLOR_BUFFER_BIT);

  if (m_simulated != m_strand.activelySimulated()) {
    m_simulated = m_strand.activelySimulated();
    m_colorChanged = true;
  }

  //    if ( m_immune != m_strand.isImmune() )
  //    {
  //        m_immune = m_strand.isImmune();
  //        m_colorChanged = true;
  //    }

  // assigns different color to clump center
  if (m_strand.getIsClumpCenterLine()) {
    m_palette["first material"] = Color(255, 0, 0);
    m_palette["second material"] = Color(0, 255, 0);
    m_palette["simple tip"] = Color(155, 200, 100);
  } else if (!m_simulated) {
    m_palette["first material"] = Color(155, 125, 75);
    m_palette["second material"] = Color(55, 125, 155);
    m_palette["simple tip"] = Color(100, 50, 75);
  } else if (m_selected) {
    m_palette["first material"] = Color(255, 175, 175);
    m_palette["second material"] = Color(155, 225, 255);
    m_palette["simple tip"] = Color(125, 195, 255);
  } else if (m_immune) {
    m_palette["first material"] = Color(255, 175, 175);
    m_palette["second material"] = Color(155, 225, 255);
    m_palette["simple tip"] = Color(125, 195, 255);
  } else  // back to original colors
  {
    m_palette["first material"] =
        Color(139, 69, 19);  // Color( 55, 125, 155 ); //Color( 255, 125, 75 );
    m_palette["second material"] =
        Color(139, 69, 19);  // Color( 55, 125, 155 ); //Color( 55, 125, 255 );
    m_palette["simple tip"] = m_palette["default simple tip"];
  }

  pushTransformationMatrix();

  if (m_drawStrand) {
    if (m_drawMode == LINES) {
      drawSimpleStrand();
    } else {
      drawSmoothStrand();
    }
  }

  if (m_drawFlow) {
    drawFlow();
  }

  if (m_drawRootMaterial) {
    drawRootMaterialFrame();
  }

  switch (m_drawForce) {
    case GRAVITATION:
      drawForce<GravitationForce>();
      break;
    case BENDING:
      drawForce<BendingForce<> >();
      break;
    case BENDING_VISCOUS:
      drawForce<BendingForce<Viscous> >();
      break;
    case STRETCHING:
      drawForce<StretchingForce<> >();
      break;
    case STRETCHING_VISCOUS:
      drawForce<StretchingForce<Viscous> >();
      break;
    case TWISTING:
      drawForce<TwistingForce<> >();
      break;
    case TWISTING_VISCOUS:
      drawForce<TwistingForce<Viscous> >();
      break;
    case NONE:
    default:
      break;
  }

  if (m_drawnVolume != NO_VOLUME) {
    drawVolume();
  }

  popTransformationMatrix();
  glPopAttrib();
}

void StrandRenderer::glVertex(int vertexIdx, bool transform) const {
  if (transform) {
    const Vec3f P = m_transformationMatrix.block<3, 3>(0, 0) *
                        Vec3f::Map(&m_currentVertices[3 * vertexIdx]) +
                    m_transformationMatrix.block<3, 1>(0, 3);
    glVertex3fv(P.data());
  } else {
    glVertex3fv(&m_currentVertices[3 * vertexIdx]);
  }
}

Vec3f StrandRenderer::vertex(int vertexIdx, bool transform) const {
  if (transform) {
    return (m_transformationMatrix.block<3, 3>(0, 0) *
                Vec3f::Map(&m_currentVertices[3 * vertexIdx]) +
            m_transformationMatrix.block<3, 1>(0, 3));
  } else {
    return Vec3f::Map(&m_currentVertices[3 * vertexIdx]);
  }
}
void StrandRenderer::drawVertices(int flag) const {
  pushTransformationMatrix();

  LockGuard lock(m_geometryMutex);

  glEnableClientState(GL_VERTEX_ARRAY);

  glVertexPointer(3, GL_FLOAT, 0, &m_currentVertices[0]);
  glDrawArrays(flag, 0, m_currentVertices.size() / 3);

  // deactivate vertex arrays after drawing
  glDisableClientState(GL_VERTEX_ARRAY);

  popTransformationMatrix();
}

void StrandRenderer::drawSimpleStrand() {
  LockGuard lock(m_geometryMutex);

  const unsigned numVertices = m_currentVertices.size() / 3;

  if (m_colorChanged || m_blendedColors.size() != m_currentVertices.size()) {
    m_colorChanged = false;

    const Color::Channel* rootColor = m_palette["simple root"].data();
    const Color::Channel* tipColor = m_palette["simple tip"].data();

    m_blendedColors.clear();
    m_blendedColors.reserve(4 * numVertices);
    for (unsigned vtx = 0; vtx < numVertices; ++vtx) {
      float t = vtx / float(numVertices - 2);

      m_blendedColors.push_back((1 - t) * rootColor[0] + t * tipColor[0]);
      m_blendedColors.push_back((1 - t) * rootColor[1] + t * tipColor[1]);
      m_blendedColors.push_back((1 - t) * rootColor[2] + t * tipColor[2]);
      m_blendedColors.push_back(1.);
    }
  }

  glLineWidth(m_lineWidth);
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);

  glVertexPointer(3, GL_FLOAT, 0, &m_currentVertices[0]);
  glColorPointer(4, GL_FLOAT, 0, &m_blendedColors[0]);

  glDrawArrays(GL_LINE_STRIP, 0, numVertices);

  //        glPointSize(4.0);
  //        glDrawArrays( GL_POINTS, 0, numVertices );
  //
  if (m_drawDebug)
    if (!m_debugVertices.empty()) {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

      int transparencyIndex = 1;
      for (auto debugPts = m_debugVertices.begin();
           debugPts != m_debugVertices.end(); ++debugPts, ++transparencyIndex) {
        // If the debugDisplaySet is empty, display all available substeps;
        // otherwise display only those that are in the set
        if (!s_debugDisplaySet.empty() &&
            s_debugDisplaySet.find(transparencyIndex) ==
                s_debugDisplaySet.end()) {
          continue;
        }

        for (unsigned vtx = 0; vtx < numVertices; ++vtx)
          m_blendedColors[4 * vtx + 3] =
              transparencyIndex / ((float)m_debugVertices.size() + 1);

        glVertexPointer(3, GL_FLOAT, 0, &(*debugPts)[0]);
        glDrawArrays(GL_LINE_STRIP, 0, debugPts->size() / 3);
      }
      glDisable(GL_BLEND);
    }
  // deactivate vertex arrays after drawing
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);

  if (m_drawDebug) {
    if (!m_hitPoints.empty()) {
      glPushAttrib(GL_CURRENT_BIT | GL_COLOR);
      glColor3f(1., 1., 0.);
      glPointSize(4);
      glBegin(GL_POINTS);
      for (auto hitPoint = m_hitPoints.begin(); hitPoint != m_hitPoints.end();
           ++hitPoint) {
        glVertex3fv(hitPoint->data());
      }
      glEnd();
      glPopAttrib();
    }
  }
}

void StrandRenderer::recomputeFlowQuads() {
  LockGuard lock(m_geometryMutex);

  if (!(m_colorChanged || m_flowChanged)) return;
  m_flowChanged = false;
  m_colorChanged = false;

  computeFlowQuads(m_flowStrand);
}

void StrandRenderer::recomputeQuads() {
  LockGuard lock(m_geometryMutex);

  if (!(m_colorChanged || m_geometryChanged)) return;
  m_geometryChanged = false;
  m_colorChanged = false;

  computeQuads(m_smoothStrand);
}

void StrandRenderer::computeFlowQuads(QuadData& quads) const {
  const int slices = 8;

  const float angleSlice = 2. * M_PI / slices;

  const unsigned numVertices = m_currentVertices.size() / 3;

  quads.m_quadVertices.resize(3 * numVertices * slices);
  quads.m_quadColors.resize(4 * numVertices * slices);
  quads.m_quadNormals.resize(quads.m_quadVertices.size());

  quads.m_quadIndices.resize(4 * (numVertices - 1) * slices);

  //        Vec4fArray slicesColors( slices );
  //        const char* c = "water";
  //
  //        for ( int k = 0; k < slices; ++k )
  //        {
  //            const Color &co = m_palette.find( c )->second;
  //            slicesColors[k][0] = co.data()[0];
  //            slicesColors[k][1] = co.data()[1];
  //            slicesColors[k][2] = co.data()[2];
  //            slicesColors[k][3] = 0.2f;
  //        }

  int num_components = m_component_colors.size();
  const VecXx& components = m_strand.getStepper()->flowComponents();

  Vec3f tangent, normal;
  unsigned curIndex = 0, curColIndex = 0;
  for (unsigned j = 0; j < numVertices; ++j) {
    // Radii

    const Vec3f& vertex = Vec3f::Map(&m_currentVertices[3 * j]);

    if (j + 1 < numVertices) {
      normal = Vec3f::Map(&m_currentMaterialFrames[3 * j]);
      tangent =
          (Vec3f::Map(&m_currentVertices[3 * (j + 1)]) - vertex).normalized();
    }

    float area = (float)(m_strand.getCurrentAreaDegreesOfFreedom()(j));
    float ra = m_strand.getParameters().getRadiusA(j, numVertices);
    float rb = m_strand.getParameters().getRadiusB(j, numVertices);
    float h = cyl_h_from_area(ra, rb, area);

    Vec3x color = Vec3x::Zero();
    for (int k = 0; k < num_components; ++k) {
      color += m_component_colors[k] * components(j * num_components + k);
    }

    for (int k = 0; k < slices; ++k) {
      // For psychedelic results: replace curColIndex by curIndex in line below
      Vec4f::Map(&quads.m_quadColors[curColIndex]) =
          Vec4f((float)color(0), (float)color(1), (float)color(2), 0.2f);

      const float angle = k * angleSlice + M_PI / 2;

      float r =
          m_strand.getParameters().getRadiusAtAngle(j, numVertices, angle);

      Vec3f::Map(&quads.m_quadNormals[curIndex]) = normal;
      Vec3f::Map(&quads.m_quadVertices[curIndex]) = vertex + (r + h) * normal;

      rotateAxisAngle(normal, tangent, angleSlice);

      curIndex += 3;
      curColIndex += 4;
    }
  }

  curIndex = 0;
  unsigned i = 0;

  for (unsigned j = 0; j + 1 < numVertices; ++j) {
    for (int k = 0; k < slices; ++k) {
      const unsigned k1 = (k + 1) % slices;
      quads.m_quadIndices[i++] = curIndex + k;
      quads.m_quadIndices[i++] = curIndex + k + slices;
      quads.m_quadIndices[i++] = curIndex + k1 + slices;
      quads.m_quadIndices[i++] = curIndex + k1;
    }

    curIndex += slices;
  }
}

void StrandRenderer::computeQuads(QuadData& quads, DrawVolume volume) const {
  float alphaCoeff = 1.f;
  switch (volume) {
    case SELF_COLLISIONS:
      alphaCoeff = m_palette.find("self collisions volume")->second.data()[3];
      break;
    case EXTERNAL_COLLISIONS:
      alphaCoeff =
          m_palette.find("external collisions volume")->second.data()[3];
      break;
    default:
      break;
  }

  const float scale = 10.0f * m_lineWidth;
  const int slices = 8;

  const float angleSlice = 2. * M_PI / slices;

  const unsigned numVertices = m_currentVertices.size() / 3;

  quads.m_quadVertices.resize(3 * numVertices * slices);
  quads.m_quadColors.resize(4 * numVertices * slices);
  quads.m_quadNormals.resize(quads.m_quadVertices.size());

  quads.m_quadIndices.resize(4 * (numVertices - 1) * slices);

  Vec4fArray slicesColors(slices);
  for (int k = 0; k < slices; ++k) {
    const char* c = NULL;

    switch (volume) {
      case SELF_COLLISIONS:
        c = "self collisions volume";
        break;
      case EXTERNAL_COLLISIONS:
        c = "external collisions volume";
        break;
      default:
        c = ((k < (slices + 3) / 4) ||
             (slices / 2 <= k && k < (3 * slices + 3) / 4))
                ? "first material"
                : "second material";
        break;
    }

    c = "random material";

    const Color& co = m_palette.find(c)->second;
    slicesColors[k][0] = co.data()[0];
    slicesColors[k][1] = co.data()[1];
    slicesColors[k][2] = co.data()[2];
    slicesColors[k][3] = alphaCoeff;
  }

  Vec3f tangent, normal;
  unsigned curIndex = 0, curColIndex = 0;
  for (unsigned j = 0; j < numVertices; ++j) {
    // Radii

    const Vec3f& vertex = Vec3f::Map(&m_currentVertices[3 * j]);

    if (j + 1 < numVertices) {
      normal = Vec3f::Map(&m_currentMaterialFrames[3 * j]);
      tangent =
          (Vec3f::Map(&m_currentVertices[3 * (j + 1)]) - vertex).normalized();
    }

    for (int k = 0; k < slices; ++k) {
      // For psychedelic results: replace curColIndex by curIndex in line below
      Vec4f::Map(&quads.m_quadColors[curColIndex]) = slicesColors[k];

      const float angle = k * angleSlice + M_PI / 2;

      //            float r ;
      //            switch ( volume )
      //            {
      //            case SELF_COLLISIONS:
      //                r = m_strand.collisionParameters().selfCollisionsRadius(
      //                j, angle ); break;
      //            case EXTERNAL_COLLISIONS:
      //                r =
      //                m_strand.collisionParameters().externalCollisionsRadius(
      //                j, angle ); break;
      //            default:
      //                r = scale * m_strand.getParameters().getRadiusAtAngle(
      //                j, angle ); break;
      //            }

      float r = m_strand.getParameters().getRadiusAtAngle(
          j, m_strand.getNumVertices(), angle);

      Vec3f::Map(&quads.m_quadNormals[curIndex]) = normal;
      Vec3f::Map(&quads.m_quadVertices[curIndex]) = vertex + r * normal;

      rotateAxisAngle(normal, tangent, angleSlice);

      curIndex += 3;
      curColIndex += 4;
    }
  }

  curIndex = 0;
  unsigned i = 0;

  for (unsigned j = 0; j + 1 < numVertices; ++j) {
    for (int k = 0; k < slices; ++k) {
      const unsigned k1 = (k + 1) % slices;
      quads.m_quadIndices[i++] = curIndex + k;
      quads.m_quadIndices[i++] = curIndex + k + slices;
      quads.m_quadIndices[i++] = curIndex + k1 + slices;
      quads.m_quadIndices[i++] = curIndex + k1;
    }

    curIndex += slices;
  }
}

void StrandRenderer::drawQuads(const QuadData& quads) {
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);

  glVertexPointer(3, GL_FLOAT, 0, &quads.m_quadVertices[0]);
  glColorPointer(4, GL_FLOAT, 0, &quads.m_quadColors[0]);
  glNormalPointer(GL_FLOAT, 0, &quads.m_quadNormals[0]);

  // TODO ? We could save some memory if we used GL_QUAD_STRIP with
  // glMultiDrawElements
  glDrawElements(GL_QUADS, quads.m_quadIndices.size(), GL_UNSIGNED_INT,
                 &quads.m_quadIndices[0]);

  // deactivate vertex arrays after drawing
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
}

void StrandRenderer::drawFlow() {
  recomputeFlowQuads();

  glEnable(GL_LIGHTING);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

  // Hack to make it work with maya's default two-sided lighting
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  drawQuads(m_flowStrand);

  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);

  glDisable(GL_BLEND);
}

void StrandRenderer::drawSmoothStrand() {
  // std::cout<<"drawing Smooth strand"<<std::endl;
  recomputeQuads();

  if (m_drawMode >= QUADS_SHADED) {
    glEnable(GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    // Hack to make it work with maya's default two-sided lighting
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

    if (m_drawMode == QUADS_BLENDED) {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
  }

  drawQuads(m_smoothStrand);

  if (m_drawMode >= QUADS_SHADED) {
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);

    if (m_drawMode == QUADS_BLENDED) {
      glDisable(GL_BLEND);
    }
  }
}

void StrandRenderer::drawVolume() {
  // std::cout<<"drawing Smooth strand"<<std::endl;

  {
    LockGuard lock(m_geometryMutex);
    computeQuads(m_volume, m_drawnVolume);
  }

  // Volume

  glEnable(GL_LIGHTING);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

  // Hack to make it work with maya's default two-sided lighting
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  drawQuads(m_volume);

  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);
  glDisable(GL_BLEND);
}

template <typename ForceT>
void StrandRenderer::drawForce() {
  static const double SCALE = 10.0;

  // FIXME m_geometryMutex only protects m_currentVertices and
  // m_currentMaterialFrames Access to other parts of the strands are *NOT*
  // thread-safe
  LockGuard lock(m_geometryMutex);
  glPushAttrib(GL_COLOR);

  VecXx force(m_strand.getNumVertices() * 4 - 1);

  // Here is a clever/ugly hack: simply doing accumulateCurrentF would give zero
  // for viscous forces because they take the current state as "rest shape".
  // Actually, if you try to call accumulateCurrentF on a dissipative force it
  // won't compile, see why in ForceAccumulator
  const_cast<ElasticStrand&>(m_strand).swapStates();
  ForceAccumulator<ForceT>::accumulateFuture(
      force, const_cast<ElasticStrand&>(m_strand));
  const_cast<ElasticStrand&>(m_strand).swapStates();

  glLineWidth(2);
  glBegin(GL_LINES);
  glColor3dv(m_palette["force"].data());

  for (int vtx = 0; vtx < m_strand.getNumVertices(); ++vtx) {
    const Eigen::Vector3d x0 = m_strand.getVertex(vtx);
    glVertex3dv(x0.data());
    const Eigen::Vector3d x1 = x0 + SCALE * force.segment<3>(4 * vtx);
    glVertex3dv(x1.data());
  }

  glEnd();
  glPopAttrib();
}

void StrandRenderer::drawRootMaterialFrame() {
  // std::cout<<"draw root material Frame"<<std::endl;
  glLineWidth(2);
  glBegin(GL_LINES);

  const Color& color1 = m_palette["first material"];
  const Color& color2 = m_palette["second material"];

  Scalar r = 1.0;

  const Vec3f& v0 = Vec3f::Map(&m_currentVertices[0]);
  const Vec3f& v1 = Vec3f::Map(&m_currentVertices[3]);
  const Vec3f& tangent = (v1 - v0).normalized();
  const Vec3f& normal = Vec3f::Map(&m_currentMaterialFrames[0]);

  const Vec3f x = (v0 + v1) * 0.5;
  //    if (m_scaleToRadius) r = m_rod.radiusA(eh);

  OpenGL::color(color1);
  const Vec3f& y1 = x + r * normal;
  glVertex3fv(x.data());
  glVertex3fv(y1.data());

  OpenGL::color(color2);

  const Vec3f& y2 = x + r * (tangent.cross(normal));
  glVertex3fv(x.data());
  glVertex3fv(y2.data());

  glEnd();
}

Vec3x StrandRenderer::calculateObjectCenter() {
  Vec3x center = Vec3x::Zero();

  const int numVertices = m_currentVertices.size() / 3;
  for (int vtx = 0; vtx < numVertices; ++vtx) {
    const Vec3x& vertex =
        Vec3f::Map(&m_currentVertices[3 * vtx]).cast<Scalar>();
    center += vertex;
  }

  center /= numVertices;

  return center;
}

Scalar StrandRenderer::calculateObjectBoundingRadius(const Vec3x& center) {
  Scalar radius = 0.0;

  const unsigned numVertices = m_currentVertices.size() / 3;
  for (unsigned vtx = 0; vtx < numVertices; ++vtx) {
    const Vec3x& vertex =
        Vec3f::Map(&m_currentVertices[3 * vtx]).cast<Scalar>();
    radius = std::max(radius, (vertex - center).norm());
  }

  return radius;
}

void StrandRenderer::pushTransformationMatrix() const {
  glPushMatrix();
  glMultMatrixf(m_transformationMatrix.data());
}

void StrandRenderer::pushInverseTransformationMatrix() const {
  glPushMatrix();
  Eigen::Matrix4f inv(m_transformationMatrix.inverse());
  glMultMatrixf(inv.data());
}

void StrandRenderer::popTransformationMatrix() const { glPopMatrix(); }

}  // namespace strandsim
