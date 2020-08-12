/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef COLLISIONBASE_HH_
#define COLLISIONBASE_HH_

#include <set>

#include "../Core/Definitions.hh"

namespace strandsim {

class ElasticStrand;

class CollisionBase {
 public:
  CollisionBase();
  virtual ~CollisionBase();

  virtual bool analyse() = 0;
  virtual bool notIn(
      const std::set<const ElasticStrand*>& alreadySeenStrands) = 0;
  virtual void putIn(std::set<const ElasticStrand*>& alreadySeenStrands) = 0;

  friend std::ostream& operator<<(std::ostream& os,
                                  const CollisionBase& collision);

 protected:
  virtual void print(std::ostream& os) const;
};

class FaceProxy;

class FaceCollision {
 public:
  virtual ~FaceCollision() {}

  virtual const FaceProxy* face() const = 0;

  Scalar faceFrictionCoefficient() const;
  virtual Scalar faceAdhesionForce() const = 0;
  virtual Scalar faceYield() const = 0;
  virtual Scalar faceEta() const = 0;
  virtual Scalar facePower() const = 0;
  virtual bool doSOCSolve() const = 0;

  Scalar m_u, m_v, m_w;
};

bool compareCT(const CollisionBase* ct1, const CollisionBase* ct2);

bool sameCT(const CollisionBase* c1, const CollisionBase* c2);

} /* namespace strandsim */
#endif /* COLLISIONBASE_HH_ */
