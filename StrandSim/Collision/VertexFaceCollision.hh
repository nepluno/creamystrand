/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef VERTEXFACECOLLISION_HH_
#define VERTEXFACECOLLISION_HH_

#include "ContinuousTimeCollision.hh"
#include "ElementProxy.hh"

namespace strandsim {

class VertexFaceCollision : public ContinuousTimeCollision,
                            public FaceCollision {
 public:
  VertexFaceCollision(ElasticStrand* firstStrand, int firstVertex,
                      const FaceProxy* faceProxy)
      : ContinuousTimeCollision(firstStrand, firstVertex),
        m_faceProxy(faceProxy),
        m_yield(0.),
        m_eta(0.),
        m_power(1.),
        m_adhesive_force(0.) {}

  virtual ~VertexFaceCollision() {}

  virtual bool analyse();

  friend bool compare(const VertexFaceCollision* vf1,
                      const VertexFaceCollision* vf2);

  Vec3x meshVelocity(Scalar dt) const { return m_meshDisplacement / dt; }

  const FaceProxy* face() const { return m_faceProxy; }

  Scalar faceAdhesionForce() const { return m_adhesive_force; }

  Scalar faceYield() const { return m_yield; }

  Scalar faceEta() const { return m_eta; }

  Scalar facePower() const { return m_power; }

  bool doSOCSolve() const { return m_do_soc_solve; }

 protected:
  void print(std::ostream& os) const;

  const FaceProxy* const m_faceProxy;
  Vec3x m_meshDisplacement;
  Vec3x m_collisionOffset;

  Scalar m_adhesive_force;
  Scalar m_yield;
  Scalar m_eta;
  Scalar m_power;
  bool m_do_soc_solve;
};

}  // namespace strandsim

#endif
