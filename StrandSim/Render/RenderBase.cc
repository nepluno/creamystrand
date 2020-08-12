/**
 * \copyright 2009 Miklos Bergou
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "RenderBase.hh"

#include "../Core/Definitions.hh"
#include "OpenGLHeaders.hh"

namespace strandsim {

void RenderBase::render() {
  if (enableNormals()) setNormalArray();
  if (enableColors()) setColorArray();
  if (enableVertices()) drawVertices();
  if (enableNormals()) unsetNormalArray();
  if (enableColors()) unsetColorArray();
}

void RenderBase::drawVertices() {
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_COLOR_MATERIAL);
  glColor3ub(255, 255, 255);

  // activate and specify pointer to vertex array
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, &(vertices[0]));

  // draw the triangles
  glDrawElements(GL_TRIANGLES, (int)vertices.size() / 3, GL_UNSIGNED_INT,
                 &(indices[0]));

  // deactivate vertex array
  glDisableClientState(GL_VERTEX_ARRAY);

  glDisable(GL_COLOR_MATERIAL);
}

void RenderBase::setNormalArray() {
  glEnableClientState(GL_NORMAL_ARRAY);
  glNormalPointer(GL_FLOAT, 0, &(normals[0]));
}

void RenderBase::unsetNormalArray() { glDisableClientState(GL_NORMAL_ARRAY); }

void RenderBase::setColorArray() {
  glEnableClientState(GL_COLOR_ARRAY);
  glColorPointer(3, GL_UNSIGNED_BYTE, 0, &(colors[0]));
}

void RenderBase::unsetColorArray() { glDisableClientState(GL_COLOR_ARRAY); }

}  // namespace strandsim
