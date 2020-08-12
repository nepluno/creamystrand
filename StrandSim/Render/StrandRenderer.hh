/**
 * \copyright 2011 Susan Howard, Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDRENDERER_HH
#define STRANDRENDERER_HH

#include <set>
#include <string>

#include "../Utils/ThreadUtils.hh"
#include "Color.hh"
#include "OpenGLDecl.hh"
#include "RenderBase.hh"

namespace strandsim {

class ElasticStrand;

/** Class that implements OpenGL rendering for rods. */
class StrandRenderer : public RenderBase {
 public:
  static std::set<int> s_debugDisplaySet;

  enum DrawMode { LINES, QUADS, QUADS_SHADED, QUADS_BLENDED };

  enum DrawForce {
    NONE,
    GRAVITATION,
    BENDING,
    BENDING_VISCOUS,
    STRETCHING,
    STRETCHING_VISCOUS,
    TWISTING,
    TWISTING_VISCOUS
  };

  enum DrawVolume { NO_VOLUME, SELF_COLLISIONS, EXTERNAL_COLLISIONS };

  StrandRenderer(const ElasticStrand &strand,
                 const std::vector<float> &currentVertices,
                 const std::vector<float> &currentMaterialFrames,
                 MutexType &geometryMutex, int num_components);

  virtual Vec3x calculateObjectCenter();
  virtual Scalar calculateObjectBoundingRadius(const Vec3x &center);

  void render();
  void drawVertices(int flag = GL_LINE_STRIP) const;

  Vec3f vertex(int vertexIdx, bool transform = true) const;
  void glVertex(int vertexIdx, bool transform = true) const;

  bool &drawStrand() { return m_drawStrand; }
  bool &drawRootMaterial() { return m_drawRootMaterial; }
  bool &drawAnchors() { return m_drawAnchors; }
  DrawMode &drawMode() { return m_drawMode; }
  DrawVolume &drawnVolume() { return m_drawnVolume; }
  void setDrawDebug(bool drawDebug) { m_drawDebug = drawDebug; }
  void setLineWidth(GLfloat lineWidth) {
    m_lineWidth = lineWidth;
    m_geometryChanged = true;
  }

  void setSelected(bool selected) {
    if (selected != m_selected) {
      m_selected = selected;
      m_colorChanged = true;
    }
  }

  void setStrandColorInSimpleMode(const Color &i_root, const Color &i_tip) {
    m_palette["simple root"] = i_root;
    m_palette["default simple tip"] = i_tip;
  }

  void ackGeometryChange() {
    m_geometryChanged = true;
    m_flowChanged = true;
  }

  const std::vector<float> &currentVertices() const {
    return m_currentVertices;
  }

  MutexType &geometryMutex() const { return m_geometryMutex; }

  const Eigen::Matrix4f &transformationMatrix() const {
    return m_transformationMatrix;
  }

  Eigen::Matrix4f &transformationMatrix() { return m_transformationMatrix; }

  void pushTransformationMatrix() const;
  void pushInverseTransformationMatrix() const;
  void popTransformationMatrix() const;

  std::list<std::vector<float> > m_debugVertices;
  std::list<Vec3f> m_hitPoints;

 protected:
  struct QuadData {
    std::vector<GLfloat> m_quadVertices;
    std::vector<GLfloat> m_quadColors;
    std::vector<GLfloat> m_quadNormals;
    std::vector<GLuint> m_quadIndices;
  };

  void drawSimpleStrand();
  void drawSmoothStrand();
  void drawFlow();
  void drawRootMaterialFrame();
  template <typename ForceT>
  void drawForce();
  void drawVolume();

  void recomputeFlowQuads();
  void recomputeQuads();
  void computeQuads(QuadData &quads, DrawVolume volume = NO_VOLUME) const;
  void computeFlowQuads(QuadData &quads) const;
  static void drawQuads(const QuadData &quads);

  const ElasticStrand &m_strand;
  const std::vector<float> &m_currentVertices;
  const std::vector<float> &m_currentMaterialFrames;

  bool m_selected;
  bool m_simulated;
  bool m_immune;
  bool m_drawStrand;
  bool m_drawFlow;
  bool m_drawRootMaterial;
  bool m_drawAnchors;
  DrawMode m_drawMode;
  DrawForce m_drawForce;
  DrawVolume m_drawnVolume;

  //    std::vector<Color> m_simpleStrand;

  std::map<std::string, Color> m_palette;

  std::vector<GLfloat> m_blendedColors;
  QuadData m_smoothStrand;
  QuadData m_flowStrand;
  QuadData m_volume;

  volatile bool m_geometryChanged;
  volatile bool m_flowChanged;
  bool m_colorChanged;
  MutexType &m_geometryMutex;
  bool m_drawDebug;
  GLfloat m_lineWidth;

  Eigen::Matrix4f m_transformationMatrix;

  std::vector<Vec3x> m_component_colors;
};

}  // namespace strandsim

#endif  // RODRENDERER_HH
