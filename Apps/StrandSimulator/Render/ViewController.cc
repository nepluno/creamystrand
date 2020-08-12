/**
 * \copyright 2012 Adrian Secord, Miklous Bergou
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "ViewController.hh"

const Scalar ViewController::eps = std::numeric_limits<Scalar>::epsilon();

ViewController::ViewController()
    : m_centerMode(CENTER_OBJECT),
      m_trackball(&m_camera),
      m_translator(&m_camera),
      m_zoomer(&m_camera) {
  m_camera.setPerspective(60, 1);
  m_camera.setViewport(100, 100);
}

void ViewController::ApplyCamera() {
  m_camera.applyViewport();
  m_camera.applyProjection();

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  m_camera.applyCamera();
}

/*
  void ViewController::setShading() const {
  static const bool extraLights = true;
  if (lightingMode) {
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  if( extraLights )
  {
  glEnable(GL_LIGHT1);

  GLfloat four[] = { 9, 9, .5, 1 };
  glLightfv(GL_LIGHT1, GL_POSITION, four);

  //four[0] = four[1] = four[2] = four[3] = 1.;
  four[0] = four[1] = four[2] = four[3] = 0.2;    // Not so bright
  glLightfv(GL_LIGHT1, GL_DIFFUSE, four);
  glLightfv(GL_LIGHT1, GL_SPECULAR, four);
  //glLightf(GL_LIGHT1, GL_SPOT_EXPONENT, 5.);
  glLightf(GL_LIGHT1, GL_SPOT_EXPONENT, 25.);
  }

  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

  } else {
  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);
  if( extraLights )
  {
  glDisable(GL_LIGHT1);
  }
  glDisable(GL_COLOR_MATERIAL);
  }

  switch (shadingMode) {
  case SHADE_FLAT:
  glShadeModel(GL_FLAT);
  break;

  case SHADE_SMOOTH:
  glShadeModel(GL_SMOOTH);
  break;

  default:
  assert(!"Unknown shading mode");
  break;
  }

  }
*/
void ViewController::setCamera(const Camera& c) { m_camera = c; }

const Camera& ViewController::getCamera() const { return m_camera; }

Camera& ViewController::getCamera() { return m_camera; }

void ViewController::setCenterMode(const CenterMode m) {
  m_centerMode = m;
  switch (m_centerMode) {
    case CENTER_WORLD_ORIGIN:
      m_camera.setViewCenter(Vec3d(0, 0, 0));
      break;

    case CENTER_OBJECT:
      m_camera.setViewCenter(m_objCenter);
      break;
  }
}

void ViewController::setViewCenter(const Vec3d& p) {
  m_objCenter = p;
  m_camera.setViewCenter(p);
}

void ViewController::getCenterOfInterest(ViewController::CenterMode* m,
                                         Vec3d* p) const {
  if (m) *m = m_centerMode;

  if (p) *p = m_camera.getViewCenter();
}

void ViewController::setViewDirection(const Vec3d& d) {
  m_camera.setEye(m_camera.getViewCenter() - d);
  m_camera.setUp(Vec3d(0, 1, 0));
}

void ViewController::setBoundingRadius(Scalar r) {
  m_boundingRadius = r;
  m_camera.setDefault3D(m_boundingRadius);
  m_translator.setScale(2 * m_boundingRadius);
  m_zoomer.setScale(2 * m_boundingRadius);
}

Scalar ViewController::getBoundingRadius() { return m_boundingRadius; }

void ViewController::beginRotationDrag(const Scalar x, const Scalar y) {
  m_trackball.start(Vec2d(x, y));
}

void ViewController::endRotationDrag(const Scalar x, const Scalar y) {
  m_trackball.update(Vec2d(x, y));
  m_trackball.stop();
}

void ViewController::beginTranslationDrag(const Scalar x, const Scalar y) {
  m_translator.start(Vec2d(x, y));
}

void ViewController::endTranslationDrag(const Scalar x, const Scalar y) {
  m_translator.update(Vec2d(x, y));
  m_translator.stop();
}

void ViewController::beginZoomDrag(const Scalar x, const Scalar y) {
  m_zoomer.start(Vec2d(x, y));
}

void ViewController::endZoomDrag(const Scalar x, const Scalar y) {
  m_zoomer.update(Vec2d(x, y));
  m_zoomer.stop();
}

void ViewController::updateDrag(const Scalar x, const Scalar y) {
  const Vec2d v(x, y);
  m_trackball.update(v);
  m_translator.update(v);
  m_zoomer.update(v);
}

void ViewController::calcTranslation(const Vec2d& start, const Vec2d& stop,
                                     const Vec3d& viewCenter, const Vec3d& eye,
                                     const Vec3d& up, const Scalar scale,
                                     Vec3d& translation) {
  // Points to the object
  const Vec3d n = (viewCenter - eye).normalized();
  // Spans the view plane in world coords
  const Vec3d w2 = (up - up.dot(n) * n).normalized();
  const Vec3d w1 = n.cross(w2);

  const Vec2d v = (stop - start).normalized() * scale;

  translation = v[0] * w1 + v[1] * w2;

  assert(translation.dot(n) < 10 * eps);
}
