/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef COLLISIONPARAMETERS_HH_
#define COLLISIONPARAMETERS_HH_

#ifndef M_PI
#define M_PI 3.141592653589793238462650288
#endif

#include "../Dynamic/CohesionTable.hh"
#include "ElasticStrandParameters.hh"

namespace strandsim {

class CollisionParameters {
 public:
  enum CollisionType { SELF, EXTERNAL };

  // For now nothing depends on edges, but maybe it will later

  CollisionParameters()
      : m_externalCollisionsRadius(.1),
        m_selfCollisionsRadius(.1),
        m_frictionCoefficient(.2),
        m_meshFrictionCoefficient(0.),
        m_ellipticalRadiiRatio(1.),
        m_cohesionSigma(VecXx::Constant(1, 72.0)),
        m_cohesionTheta(40.8 / 180.0 * M_PI),
        m_cohesionMaxDist(0.25),
        m_reactsToSelfCollisions(true),
        m_impulseMaxNorm(10000.0),
        m_createsSelfCollisions(true),
        m_fakeLayering(false),
        m_constantCollisionRadius(false),
        m_associatedStrandParams(NULL),
        m_maxNumCollisionsPerEdge(6) {}

  CollisionParameters(
      Scalar selfCollisionsRadius, Scalar externalCollisionsRadius,
      Scalar frictionCoefficient, Scalar meshFrictionCoefficient,
      Scalar ellipticalRadiiRatio, const VecXx& cohesionSigma,
      Scalar cohesionTheta, Scalar cohesionMaxDist, Scalar impulseMaxNorm,
      int maxNumCollsionsPerEdge, bool reactsToSelfCollisions,
      bool createsSelfCollisions = true, bool ellipticalCollisions = false,
      bool constantCollisionRadius = false)
      : m_externalCollisionsRadius(externalCollisionsRadius),
        m_selfCollisionsRadius(selfCollisionsRadius),
        m_frictionCoefficient(frictionCoefficient),
        m_meshFrictionCoefficient(meshFrictionCoefficient),
        m_ellipticalRadiiRatio(ellipticalRadiiRatio),
        m_reactsToSelfCollisions(reactsToSelfCollisions),
        m_createsSelfCollisions(createsSelfCollisions),
        m_cohesionSigma(cohesionSigma),
        m_cohesionTheta(cohesionTheta),
        m_cohesionMaxDist(cohesionMaxDist),
        m_impulseMaxNorm(impulseMaxNorm),
        m_maxNumCollisionsPerEdge(maxNumCollsionsPerEdge),
        m_constantCollisionRadius(constantCollisionRadius),
        m_fakeLayering(false),
        m_associatedStrandParams(NULL) {}

  void initializeCohesionTable() {
    m_cohesionTable.setParameter(1.0, m_cohesionTheta, m_selfCollisionsRadius,
                                 m_cohesionMaxDist, 1024);
    m_cohesionTable.construct_alpha_table();

    m_cohesionTable_planar.setParameter(1.0, m_cohesionTheta,
                                        m_selfCollisionsRadius,
                                        m_cohesionMaxDist * 0.5, 512);
    m_cohesionTable_planar.construct_planar_alpha_table();
  }

  bool usesFakeLayering() const { return m_fakeLayering; }

  Scalar externalCollisionsRadius(unsigned edgeIdx, Scalar angle = 0.) const {
    return collisionsRadius(EXTERNAL, edgeIdx, angle);
  }
  Scalar& externalCollisionsRadius() { return m_externalCollisionsRadius; }

  Scalar selfCollisionsRadius(unsigned edgeIdx, Scalar angle = 0.) const {
    return collisionsRadius(SELF, edgeIdx, angle);
  }
  Scalar& selfCollisionsRadius() { return m_selfCollisionsRadius; }

  Scalar collisionsRadius(CollisionType type, unsigned edgeIdx,
                          unsigned numVertices, Scalar angle = 0.) const {
    const Scalar baseRadius =
        type == SELF ? m_selfCollisionsRadius : m_externalCollisionsRadius;

    if (m_associatedStrandParams && !m_constantCollisionRadius) {
      const Scalar majorRad =
          m_associatedStrandParams->interpolatedRadiusMultiplier(
              (int)edgeIdx, (int)numVertices) *
          baseRadius / m_associatedStrandParams->getMaxRadiusMultiplier();

      return majorRad;
    }
    return baseRadius;
  }

  Scalar ellipticalRadiiRatio(unsigned /*edgeIdx*/) const {
    return m_ellipticalRadiiRatio;
  }

  Scalar frictionCoefficient(unsigned /*edgeIdx*/) const {
    return m_frictionCoefficient;
  }

  bool reactsToSelfCollisions() const { return m_reactsToSelfCollisions; }

  bool createsSelfCollisions() const { return m_createsSelfCollisions; }

  void multiplyExternalCollisionsRadius(Scalar m) {
    this->m_externalCollisionsRadius *= m;
  }

  void multiplyFrictionCoefficient(Scalar m) {
    this->m_frictionCoefficient *= m;
  }

  void multiplySelfCollisionsRadius(Scalar m) {
    this->m_selfCollisionsRadius *= m;
  }

  void useFakeLayering(bool use) { m_fakeLayering = use; }

  void setAssociatedStrandParameters(
      const ElasticStrandParameters& strandParams) {
    m_associatedStrandParams = &strandParams;
  }

  Scalar adhesionForce(const Scalar& A, const Scalar& d0,
                       const VecXx& color) const {
    return m_cohesionTable.interpolate_dEdd(A, d0) * m_cohesionSigma.dot(color);
  }

  Scalar adhesionForcePlanar(const Scalar& A, const Scalar& d0,
                             const VecXx& color) const {
    return m_cohesionTable_planar.interpolate_dEdd_planar(A, d0) *
           m_cohesionSigma.dot(color);
  }

  // private:
  Scalar m_externalCollisionsRadius;
  Scalar m_selfCollisionsRadius;  // Made public for sheer lazines
  Scalar m_frictionCoefficient;
  Scalar m_meshFrictionCoefficient;
  Scalar m_ellipticalRadiiRatio;

  VecXx m_cohesionSigma;
  Scalar m_cohesionTheta;
  Scalar m_cohesionMaxDist;
  Scalar m_impulseMaxNorm;

  int m_maxNumCollisionsPerEdge;

  bool m_reactsToSelfCollisions;
  bool m_createsSelfCollisions;
  bool m_constantCollisionRadius;
  bool m_fakeLayering;

  CohesionTable m_cohesionTable;
  CohesionTable m_cohesionTable_planar;

  // Required for getting radius interpolation
  const ElasticStrandParameters* m_associatedStrandParams;
};

}  // namespace strandsim

#endif /* COLLISIONPARAMETERS_HH_ */
