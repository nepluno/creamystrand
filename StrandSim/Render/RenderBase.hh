/**
 * \copyright 2009 Miklos Bergou
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef S_RENDERBASE_HH
#define S_RENDERBASE_HH

#include <vector>

#include "../Core/Definitions.hh"
#include "OpenGLHeaders.hh"

namespace strandsim {

/** Interface for all renderers. */
class RenderBase {
 public:
  RenderBase()
      : m_enableVertices(false), m_enableNormals(false), m_enableColors(false) {
    m_enableVertices = true;

    vertices.push_back(0);
    vertices.push_back(0);
    vertices.push_back(0);

    vertices.push_back(1);
    vertices.push_back(0);
    vertices.push_back(0);

    vertices.push_back(0.5);
    vertices.push_back(1);
    vertices.push_back(0);

    indices.push_back(0);
    indices.push_back(1);
    indices.push_back(2);
  }

  virtual ~RenderBase() {}

  virtual void render();

  bool& enableVertices() { return m_enableVertices; }
  bool& enableNormals() { return m_enableNormals; }
  bool& enableColors() { return m_enableColors; }

  // virtual Vec3d calculateObjectCenter() = 0;
  // virtual Scalar calculateObjectBoundingRadius(const Vec3d& center) = 0;

 protected:
  void setNormalArray();
  void unsetNormalArray();
  void setColorArray();
  void unsetColorArray();
  void drawVertices();

  std::vector<GLfloat> vertices;
  std::vector<GLfloat> normals;
  std::vector<GLubyte> colors;
  std::vector<GLuint> indices;

  bool m_enableVertices;
  bool m_enableNormals;
  bool m_enableColors;
};

}  // namespace strandsim

#endif  // S_RENDERBASE_HH
