/**
 * \copyright 2012 Adrian Secord, Miklous Bergou
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ZOOMER_H
#define ZOOMER_H

//#include "../Core/EigenIncludes.hh"
//#include "../Core/STLIncludes.hh"
#include "../../StrandSim/Render/OpenGLHeaders.hh"
#include "Camera.hh"

/** Converts mouse motions into world translations. Modifies the
    camera's eye and view center directly. */
class Zoomer
{
public:

  /// Default constructor
  explicit Zoomer(Camera* c, const Scalar scale = 1.0);

  /// Set a particular camera to use.
  void setCamera(Camera* c);

  /// Set the scaling for mouse motions.
  void setScale(const Scalar s);

  /// Start a mouse motion.
  /// Position in [-1,1] x [-1,1].
  void start(const Vec2d& p);

  /// Update a mouse motion with a new mouse point.
  /// Position in [-1,1] x [-1,1].
  void update(const Vec2d& p);

  /// Stop a mouse motion.
  void stop();

protected:

  Camera* m_camera;                       ///< The current camera to obey.
  bool m_translating;                     ///< Whether we are translating.
  Vec2d m_startPos;                       ///< Start of a translation drag.
  Scalar m_scale;                         ///< Scaling of mouse motions.
};


#endif // ZOOMER_H
