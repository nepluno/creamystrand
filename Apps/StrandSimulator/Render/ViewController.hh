/**
 * \copyright 2012 Adrian Secord, Miklous Bergou
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef VIEWCONTROLLER_HH
#define VIEWCONTROLLER_HH

//#include "../Core/EigenIncludes.hh"
//#include "../Core/STLIncludes.hh"
#include "../../StrandSim/Render/OpenGLHeaders.hh"
#include "Camera.hh"
#include "TrackBall.hh"
#include "Translator.hh"
#include "Zoomer.hh"

class ViewController
{
public:

  /// Center of interest options
  enum CenterMode {
    CENTER_WORLD_ORIGIN,                   ///< Center the view on the origin
    CENTER_OBJECT,                         ///< Center the view on the mesh
  };

  /// Default constructor.
  ViewController();

  void ApplyCamera();

  /// Set the current camera.
  void setCamera(const Camera& c);

  /// Get the current camera.
  const Camera& getCamera() const;

  /// Get the current camera.
  Camera& getCamera();

  /// Set the view center of interest
  void setCenterMode(const CenterMode m);
  void setViewCenter(const Vec3d& p);

  void setBoundingRadius(Scalar r);
  Scalar getBoundingRadius();

  /// Get the view center of interest.  Either parameter can be NULL if
  /// it is not needed.
  void getCenterOfInterest(CenterMode* m, Vec3d* p = NULL) const;

  /// Set a particular viewing direction
  void setViewDirection(const Vec3d& d);

  /** \name Interaction commands

      Commands to modify the camera from user interaction. */

  //@{

  /// Begin a rotation drag at a point in pixel coordinates.
  void beginRotationDrag(const Scalar x, const Scalar y);

  /// End a rotation drag at a point in pixel coordinates.
  void endRotationDrag(const Scalar x, const Scalar y);

  /// Begin a translation drag at a point in pixel coordinates.
  void beginTranslationDrag(const Scalar x, const Scalar y);

  /// End a translation drag at a point in pixel coordinates.
  void endTranslationDrag(const Scalar x, const Scalar y);

  /// Begin a zoom drag at a point in pixel coordinates.
  void beginZoomDrag(const Scalar x, const Scalar y);

  /// End a zoom drag at a point in pixel coordinates.
  void endZoomDrag(const Scalar x, const Scalar y);

  /// Update a current drag operation with new coordinates.
  void updateDrag(const Scalar x, const Scalar y);

  //@}

private:

  /// Smallest amount that can be added to a scalar and cause a difference.
  static const Scalar eps;

  /// Compute the translation needed to move a point from world-space start to
  /// intersect the perpendicular line that pierces the screen at (x,y).
  static void calcTranslation(const Vec2d& start,
                              const Vec2d& stop,
                              const Vec3d& viewCenter,
                              const Vec3d& eye,
                              const Vec3d& up,
                              const Scalar scale,
                              Vec3d& translation);

  /// Not copyable
  ViewController(const ViewController& m);

  /// Not copyable
  ViewController& operator= (const ViewController& m);

  Camera m_camera;                       ///< Camera for viewing scene.

  CenterMode m_centerMode;               ///< Centering mode

  Vec3d m_objCenter;                     ///< Center of the object.
  Scalar m_boundingRadius;               ///< Bounding radius of mesh.

  TrackBall m_trackball;                 ///< Trackball for user rotation.
  Translator m_translator;               ///< User translations.
  Zoomer m_zoomer;                       ///< User zooms.

};

#endif // VIEWCONTROLLER_HH
