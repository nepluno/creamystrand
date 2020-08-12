/*
 * (c) Copyright 1993, 1994, Silicon Graphics, Inc.
 * ALL RIGHTS RESERVED
 * Permission to use, copy, modify, and distribute this software for
 * any purpose and without fee is hereby granted, provided that the above
 * copyright notice appear in all copies and that both the copyright notice
 * and this permission notice appear in supporting documentation, and that
 * the name of Silicon Graphics, Inc. not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission.
 *
 * THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
 * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
 * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
 * GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
 * SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
 * KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
 * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
 * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
 * POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * US Government Users Restricted Rights
 * Use, duplication, or disclosure by the Government is subject to
 * restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
 * (c)(1)(ii) of the Rights in Technical Data and Computer Software
 * clause at DFARS 252.227-7013 and/or in similar or successor
 * clauses in the FAR or the DOD or NASA FAR Supplement.
 * Unpublished-- rights reserved under the copyright laws of the
 * United States.  Contractor/manufacturer is Silicon Graphics,
 * Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
 *
 * OpenGL(TM) is a trademark of Silicon Graphics, Inc.
 */
/*
 * trackball.h
 * A virtual trackball implementation
 * Written by Gavin Bell for Silicon Graphics, November 1988.
 */

#include "TrackBall.hh"

namespace {

inline void rotate(Scalar m[4][4], const Vec3d& axis, const Scalar radians) {
  Vec3d naxis(axis);
  naxis.normalize();

  Scalar x = naxis[0];
  Scalar y = naxis[1];
  Scalar z = naxis[2];

  Scalar c = cos(radians);
  Scalar s = sin(radians);
  Scalar t = 1 - c;

  m[0][0] = t * x * x + c;
  m[0][1] = t * x * y - z * s;
  m[0][2] = t * x * z + y * s;
  m[1][0] = t * x * y + z * s;
  m[1][1] = t * y * y + c;
  m[1][2] = t * y * z - x * s;
  m[2][0] = t * x * z - y * s;
  m[2][1] = t * y * z + x * s;
  m[2][2] = t * z * z + c;

  m[0][3] = 0;
  m[1][3] = 0;
  m[2][3] = 0;
  m[3][0] = 0;
  m[3][1] = 0;
  m[3][2] = 0;
  m[3][3] = 1;
}

}  // namespace

TrackBall::TrackBall(Camera* c)
    : m_camera(c), m_rotating(false), m_mode(GIMBAL) {}

void TrackBall::start(const Vec2d& p) {
  assert(p[0] >= -1 && p[0] <= 1 && p[1] >= -1 && p[1] <= 1);
  m_rotating = true;
  m_startPos = p;
  for (int i = 0; i < 4; ++i) m_rotation[i] = 0;
}

void TrackBall::update(const Vec2d& p) {
  if (!m_rotating) return;

  Scalar m[4][4];
  if (m_startPos != p) {
    const Scalar coef(M_PI / 2.0);

    const Scalar left_right_motion = p[0] - m_startPos[0];
    const Scalar up_down_motion = p[1] - m_startPos[1];

    // rotate about the 'up' vector
    Vec3d up;
    if (m_mode == GIMBAL)
      up = Vec3d(0, 1, 0);
    else
      up = m_camera->getUp();
    int sign = (up.dot(m_camera->getUp()) > 0) ? 1 : -1;
    rotate(m, up, coef * sign * left_right_motion);
    m_camera->orbit(m);

    // rotate about the horizontal vector
    Vec3d horizontal =
        m_camera->getUp().cross(m_camera->getEye() - m_camera->getViewCenter());
    rotate(m, horizontal, -coef * up_down_motion);
    m_camera->orbit(m);

    m_startPos = p;
  }
}

void TrackBall::stop() { m_rotating = false; }

void TrackBall::getRotation(Scalar r[4]) const {
  for (int i = 0; i < 4; ++i) {
    r[i] = m_rotation[i];
  }
}

TrackBall::TrackBallMode TrackBall::getTrackBallMode() const { return m_mode; }

void TrackBall::setTrackBallMode(const TrackBall::TrackBallMode m) {
  m_mode = m;
}
