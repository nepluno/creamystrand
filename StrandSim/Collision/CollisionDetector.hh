/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef COLLISIONDETECTOR_HH_
#define COLLISIONDETECTOR_HH_

#include <list>
#include <vector>

#include "../Utils/SpatialHashMapFwd.hh"
#include "BVH.hh"

namespace strandsim {

class ElementProxy;
class EdgeProxy;
class FaceProxy;
class CollisionBase;
class ElementProxySortedAABBFunctor;

class CollisionDetector {
 public:
  CollisionDetector(std::vector<ElementProxy*>& elements);
  virtual ~CollisionDetector();

  void buildBVH(bool statique);
  void updateBoundingBoxes();
  void updateBoundingBox(BVHNodeType& node);

  void findCollisions(bool ignoreContinuousTime = false,
                      bool ignoreProximity = false);

  void clear();
  bool empty();

  const std::list<CollisionBase*>& getContinuousTimeCollisions() const {
    return m_continuousTimeCollisions;
  }
  std::list<CollisionBase*>& getContinuousTimeCollisions() {
    return m_continuousTimeCollisions;
  }
  const std::list<CollisionBase*>& getProximityCollisions() const {
    return m_proximityCollisions;
  }
  std::list<CollisionBase*>& getProximityCollisions() {
    return m_proximityCollisions;
  }

  static void setMaxSizeForElementBBox(double s) {
    s_maxSizeForElementBBox = s;
  }

 protected:
  void computeCollisions(const BVHNodeType& node_a, const BVHNodeType& node_b);

  template <typename CollisionT>
  bool appendCollision(ElementProxy* elem_a, ElementProxy* elem_b);

  bool appendCollision(EdgeProxy* edge_a, EdgeProxy* edge_b);
  bool appendCollision(EdgeProxy* edge_a, const FaceProxy* triangle_b);
  bool appendIntersection(EdgeProxy* edge_a, const FaceProxy* triangle_b);

  void filterWithSpatialHashMap(const Scalar largestBBox);

  std::vector<ElementProxy*> m_elementProxies;
  BVH m_bvh;
  std::list<CollisionBase*> m_continuousTimeCollisions;
  std::list<CollisionBase*> m_proximityCollisions;

  int m_broadPhaseHitCounter;
  bool m_ignoreContinuousTime;
  bool m_ignoreProximity;
  static Scalar s_maxSizeForElementBBox;

  ElementProxySortedAABBFunctor* m_sortedAABBFunctor;
  typedef SpatialHashMap<const ElementProxySortedAABBFunctor, unsigned, true>
      SpatialHashMapT;
  SpatialHashMapT* m_hashMap;
};

} /* namespace strandsim */
#endif /* COLLISIONDETECTOR_HH_ */
