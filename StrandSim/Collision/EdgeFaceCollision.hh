/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EDGEFACECOLLISION_HH_
#define EDGEFACECOLLISION_HH_

#include "ContinuousTimeCollision.hh"
#include "ElementProxy.hh"

namespace strandsim {

class EdgeFaceCollision : public EdgeCollision, public FaceCollision {
 public:
  EdgeFaceCollision(ElasticStrand* firstStrand, int firstVertex,
                    const FaceProxy* face, int firstIdx, int secondIdx,
                    int firstApex, int secondApex)
      : EdgeCollision(firstStrand, firstVertex),  //
        m_mesh(face->getMesh()),                  //
        m_face(face),                             //
        m_firstIdx(firstIdx),                     //
        m_secondIdx(secondIdx),
        m_firstApex(firstApex),
        m_secondApex(secondApex),
        m_onBoundary(face->getFace().edgeOnBoundary(m_firstApex)),
        m_adhesive_force(0.),
        m_yield(0.),
        m_eta(0.),
        m_power(1.) {
    // Maintain canonical order for sorting
    if (m_secondIdx < m_firstIdx) {
      std::swap(m_secondIdx, m_firstIdx);
      std::swap(m_firstApex, m_secondApex);
    }
  }

  virtual bool analyse();

  friend bool compare(const EdgeFaceCollision* ef1,
                      const EdgeFaceCollision* ef2);

  Vec3x meshVelocity(Scalar dt) const { return m_meshDisplacement / dt; }

  uint32_t faceEdgeId() const {
    return m_face->uniqueId() + (4 - m_firstApex - m_secondApex);
  }

  const FaceProxy* face() const { return m_face; }

  Scalar faceAdhesionForce() const { return m_adhesive_force; }

  Scalar faceYield() const { return m_yield; }

  Scalar faceEta() const { return m_eta; }

  Scalar facePower() const { return m_power; }

  bool doSOCSolve() const { return m_do_soc_solve; }

 protected:
  void print(std::ostream& os) const;

  const std::shared_ptr<TriangularMesh> m_mesh;
  const FaceProxy* m_face;
  int m_firstIdx, m_secondIdx, m_firstApex, m_secondApex;
  bool m_onBoundary;
  Vec3x m_meshDisplacement;
  Scalar m_adhesive_force;
  Scalar m_yield;
  Scalar m_eta;
  Scalar m_power;
  bool m_do_soc_solve;
};

}  // namespace strandsim

#endif
