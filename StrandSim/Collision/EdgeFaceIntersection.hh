/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EDGEFACEINTERSECTION_HH_
#define EDGEFACEINTERSECTION_HH_

#include "CollisionBase.hh"
#include "ElementProxy.hh"

namespace strandsim {

class ElasticStrand;

class EdgeFaceIntersection : public CollisionBase, public FaceCollision {
 public:
  static double s_doProximityDetection;

  EdgeFaceIntersection(const ElasticStrand* strand, int edge, EdgeProxy* edgeP,
                       const strandsim::FaceProxy* triangle);
  virtual ~EdgeFaceIntersection();

  virtual bool analyse();

  virtual bool notIn(const std::set<const ElasticStrand*>& alreadySeenStrands) {
    return find(alreadySeenStrands.begin(), alreadySeenStrands.end(),
                m_strand) == alreadySeenStrands.end();
  }
  virtual void putIn(std::set<const ElasticStrand*>& alreadySeenStrands) {
    alreadySeenStrands.insert(m_strand);
  }

  friend std::ostream& operator<<(std::ostream& os,
                                  const CollisionBase& collision);

  EdgeProxy* getEdge() { return m_edgeP; }
  int getEdgeNumber() const { return m_edge; }
  const ElasticStrand* getStrand() const { return m_strand; }
  Scalar getEdgeAbscissa() const { return m_s; }
  const FaceProxy* face() const { return m_triangle; }

  Scalar faceAdhesionForce() const { return m_adhesive_force; }

  Scalar faceYield() const { return m_yield; }

  Scalar faceEta() const { return m_eta; }

  Scalar facePower() const { return m_power; }

  bool doSOCSolve() const { return true; }

  Vec3x getFaceNormal() const;
  Vec3x getFaceVelocity(const Scalar dt) const;

  uint32_t getFaceId() const;

 protected:
  virtual void print(std::ostream& os) const;

  bool analyseProximity();

  const ElasticStrand* const m_strand;
  const int m_edge;
  EdgeProxy* const m_edgeP;
  const FaceProxy* const m_triangle;
  Scalar m_s /*, u, v, w*/;
  Scalar m_distance;
  Scalar m_adhesive_force;
  Scalar m_yield;
  Scalar m_eta;
  Scalar m_power;
  bool m_isBoundaryedgeCollision;

  Vec3x m_normal;

 private:
  bool barycentricCoordinatesIfCloser(
      const Scalar radius, const Scalar A0, const Scalar A1, const Scalar theta,
      const Scalar s0, const Vec3x& p0, const Vec3x& p1, const Vec3x& q0,
      const Vec3x& q1, const Vec3x& q2, Scalar& dist, Scalar& s, Scalar& u,
      Scalar& v, Scalar& w, Scalar& al, Scalar& A_adh, Scalar& d_adh);
};

} /* namespace strandsim */
#endif /* EDGEFACEINTERSECTION_HH_ */
