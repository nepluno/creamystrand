/**
 * \copyright 2012 Adrian Secord, Miklous Bergou
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Zoomer.hh"

Zoomer::Zoomer(Camera* c, const Scalar s)
    : m_camera(c), m_translating(false), m_scale(s) {}

void Zoomer::setCamera(Camera* c) { m_camera = c; }

void Zoomer::setScale(const Scalar s) { m_scale = s; }

void Zoomer::start(const Vec2d& p) {
  m_translating = true;
  m_startPos = p;
}

void Zoomer::update(const Vec2d& p) {
  if (!m_translating) return;

  assert(m_camera);
  Vec3d in;
  m_camera->getSpanningSet(NULL, NULL, &in);

  const Vec3d translation = in * m_scale * (p[1] - m_startPos[1]);

  m_camera->translateEye(translation);

  m_startPos = p;
}

void Zoomer::stop() { m_translating = false; }
