/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "CollisionDetector.hh"

#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Utils/SpatialHashMap.hh"
#include "../Utils/TextLog.hh"
#include "../Utils/TimeUtils.hh"
#include "ContinuousTimeCollision.hh"
#include "EdgeFaceCollision.hh"
#include "EdgeFaceIntersection.hh"
#include "ElementProxy.hh"
#include "VertexFaceCollision.hh"

namespace strandsim {
Scalar CollisionDetector::s_maxSizeForElementBBox = 1e2;

CollisionDetector::CollisionDetector(std::vector<ElementProxy*>& elements)
    : m_elementProxies(elements),
      m_bvh(),
      m_broadPhaseHitCounter(0),
      m_ignoreContinuousTime(false),
      m_ignoreProximity(false),
      m_sortedAABBFunctor(NULL),
      m_hashMap(NULL) {}

CollisionDetector::~CollisionDetector() {
  delete m_sortedAABBFunctor;
  delete m_hashMap;
}

void CollisionDetector::buildBVH(bool statique) {
  Scalar largestElemBBoxSize = 0;

#pragma omp parallel for
  for (int elemId = 0; elemId < m_elementProxies.size(); ++elemId) {
    ElementProxy* elem = m_elementProxies[elemId];

    elem->updateBoundingBox(statique);

    const BBoxType& elemBBox = elem->getBoundingBox();
    const Scalar elemBBoxSize = elemBBox.maxDim();

    if (elemBBoxSize > s_maxSizeForElementBBox) {
      WarningStream(g_log, "")
          << "Element " << *elem << " has large bounding box: " << elemBBox
          << "; will not be considered for collision detection";
      elem->resetBoundingBox();
    }

#pragma omp flush(largestElemBBoxSize)
    if (elemBBoxSize > largestElemBBoxSize) {
#pragma omp critical
      {
        if (elemBBoxSize > largestElemBBoxSize)
          largestElemBBoxSize = elemBBoxSize;
      }
    }
  }

  DebugStream(g_log, "")
      << "CollisionDetector::buildBVH(): largest element BBox size = "
      << largestElemBBoxSize;

  filterWithSpatialHashMap(largestElemBBoxSize);

  ElementProxyBBoxFunctor bboxfunctor(m_elementProxies);
  BVHBuilder<ElementProxyBBoxFunctor> bvh_builder;
  bvh_builder.build(bboxfunctor, &m_bvh);
}

void CollisionDetector::updateBoundingBoxes() {
  updateBoundingBox(m_bvh.GetNode(0));
}

void CollisionDetector::updateBoundingBox(BVHNodeType& node) {
  BVHNodeType::BBoxType& bbox = node.BBox();
  bbox.reset();
  if (node.IsLeaf())  // The leaf's bounding box contains the whole trajectory
                      // of its object(s) during this time step.
  {
    const uint32_t leaf_begin = node.LeafBegin();
    const uint32_t leaf_end = node.LeafEnd();
    for (uint32_t i = leaf_begin; i < leaf_end; ++i) {
      m_elementProxies[i]->updateBoundingBox(false);
      const BBoxType& elemBBox = m_elementProxies[i]->getBoundingBox();

      if (elemBBox.maxDim() <= s_maxSizeForElementBBox) {
        bbox.insert(elemBBox);
      } else {
        WarningStream(g_log, "")
            << "Element " << *m_elementProxies[i]
            << " has large bounding box: " << elemBBox
            << "; will not be considered for collision detection";
      }
    }
  } else  // Update the children, then this node's bounding box
  {
    BVHNodeType& hansel = m_bvh.GetNode(node.ChildIndex());
    updateBoundingBox(hansel);
    bbox.insert(hansel.BBox());
    BVHNodeType& gretel = m_bvh.GetNode(node.ChildIndex() + 1);
    updateBoundingBox(gretel);
    bbox.insert(gretel.BBox());
  }
}

void CollisionDetector::findCollisions(bool ignoreContinuousTime,
                                       bool ignoreProximity) {
  assert(empty());

  m_ignoreContinuousTime = ignoreContinuousTime;
  m_ignoreProximity = ignoreProximity;

  m_broadPhaseHitCounter = 0;
  const BVHNodeType& root = m_bvh.GetNode(0);

  if (root.IsLeaf())  // Can't really call this a tree, can we?
  {
    computeCollisions(root, root);
    return;
  }

  const BVHNodeType& h = m_bvh.GetNode(root.ChildIndex());
  const BVHNodeType& g = m_bvh.GetNode(root.ChildIndex() + 1);

  // If tree has depth 1, detect collisions at this level.
  if (h.IsLeaf() || g.IsLeaf()) {
#pragma omp parallel sections
    {
#pragma omp section
      { computeCollisions(h, h); }
#pragma omp section
      { computeCollisions(h, g); }
#pragma omp section
      { computeCollisions(g, g); }
    }
    return;
  }

  const BVHNodeType& hh = m_bvh.GetNode(h.ChildIndex());
  const BVHNodeType& hg = m_bvh.GetNode(h.ChildIndex() + 1);
  const BVHNodeType& gh = m_bvh.GetNode(g.ChildIndex());
  const BVHNodeType& gg = m_bvh.GetNode(g.ChildIndex() + 1);

  // If tree has depth 2, detect collisions at this level.
  if (hh.IsLeaf() || hg.IsLeaf() || gh.IsLeaf() || gg.IsLeaf()) {
#pragma omp parallel sections
    {
#pragma omp section
      { computeCollisions(hh, hh); }
#pragma omp section
      { computeCollisions(hh, hg); }
#pragma omp section
      { computeCollisions(hh, gh); }
#pragma omp section
      { computeCollisions(hh, gg); }
#pragma omp section
      { computeCollisions(hg, hg); }
#pragma omp section
      { computeCollisions(hg, gh); }
#pragma omp section
      { computeCollisions(hg, gg); }
#pragma omp section
      { computeCollisions(gh, gh); }
#pragma omp section
      { computeCollisions(gh, gg); }
#pragma omp section
      { computeCollisions(gg, gg); }
    }
    return;
  }

  // If the tree is deep enough, launch recursive parallel collision detection
  // from this level.
  const BVHNodeType& hhh = m_bvh.GetNode(hh.ChildIndex());
  const BVHNodeType& hhg = m_bvh.GetNode(hh.ChildIndex() + 1);
  const BVHNodeType& hgh = m_bvh.GetNode(hg.ChildIndex());
  const BVHNodeType& hgg = m_bvh.GetNode(hg.ChildIndex() + 1);
  const BVHNodeType& ghh = m_bvh.GetNode(gh.ChildIndex());
  const BVHNodeType& ghg = m_bvh.GetNode(gh.ChildIndex() + 1);
  const BVHNodeType& ggh = m_bvh.GetNode(gg.ChildIndex());
  const BVHNodeType& ggg = m_bvh.GetNode(gg.ChildIndex() + 1);

#pragma omp parallel sections
  {
#pragma omp section
    { computeCollisions(hhh, hhh); }
#pragma omp section
    { computeCollisions(hhh, hhg); }
#pragma omp section
    { computeCollisions(hhh, hgh); }
#pragma omp section
    { computeCollisions(hhh, hgg); }
#pragma omp section
    { computeCollisions(hhh, ghh); }
#pragma omp section
    { computeCollisions(hhh, ghg); }
#pragma omp section
    { computeCollisions(hhh, ggh); }
#pragma omp section
    { computeCollisions(hhh, ggg); }
#pragma omp section
    { computeCollisions(hhg, hhg); }
#pragma omp section
    { computeCollisions(hhg, hgh); }
#pragma omp section
    { computeCollisions(hhg, hgg); }
#pragma omp section
    { computeCollisions(hhg, ghh); }
#pragma omp section
    { computeCollisions(hhg, ghg); }
#pragma omp section
    { computeCollisions(hhg, ggh); }
#pragma omp section
    { computeCollisions(hhg, ggg); }
#pragma omp section
    { computeCollisions(hgh, hgh); }
#pragma omp section
    { computeCollisions(hgh, hgg); }
#pragma omp section
    { computeCollisions(hgh, ghh); }
#pragma omp section
    { computeCollisions(hgh, ghg); }
#pragma omp section
    { computeCollisions(hgh, ggh); }
#pragma omp section
    { computeCollisions(hgh, ggg); }
#pragma omp section
    { computeCollisions(hgg, hgg); }
#pragma omp section
    { computeCollisions(hgg, ghh); }
#pragma omp section
    { computeCollisions(hgg, ghg); }
#pragma omp section
    { computeCollisions(hgg, ggh); }
#pragma omp section
    { computeCollisions(hgg, ggg); }
#pragma omp section
    { computeCollisions(ghh, ghh); }
#pragma omp section
    { computeCollisions(ghh, ghg); }
#pragma omp section
    { computeCollisions(ghh, ggh); }
#pragma omp section
    { computeCollisions(ghh, ggg); }
#pragma omp section
    { computeCollisions(ghg, ghg); }
#pragma omp section
    { computeCollisions(ghg, ggh); }
#pragma omp section
    { computeCollisions(ghg, ggg); }
#pragma omp section
    { computeCollisions(ggh, ggh); }
#pragma omp section
    { computeCollisions(ggh, ggg); }
#pragma omp section
    { computeCollisions(ggg, ggg); }
  }
}

void CollisionDetector::clear() {
  for (auto collision = m_continuousTimeCollisions.begin();
       collision != m_continuousTimeCollisions.end(); ++collision) {
    delete *collision;
  }
  for (auto collision = m_proximityCollisions.begin();
       collision != m_proximityCollisions.end(); ++collision) {
    delete *collision;
  }

  m_continuousTimeCollisions.clear();
  m_proximityCollisions.clear();
}

bool CollisionDetector::empty() {
  return m_continuousTimeCollisions.empty() && m_proximityCollisions.empty();
}

template <>
bool CollisionDetector::appendCollision<ContinuousTimeCollision>(
    ElementProxy* elem_a, ElementProxy* elem_b) {
  EdgeProxy* const edge_a = dynamic_cast<EdgeProxy*>(elem_a);
  EdgeProxy* const edge_b = dynamic_cast<EdgeProxy*>(elem_b);

  if (edge_a && edge_b) {
    return false;
  } else if (edge_a) {
    FaceProxy* const triangle_b = dynamic_cast<FaceProxy*>(elem_b);
    return triangle_b && appendCollision(edge_a, triangle_b);
  } else if (edge_b) {
    FaceProxy* const triangle_a = dynamic_cast<FaceProxy*>(elem_a);
    return triangle_a && appendCollision(edge_b, triangle_a);
  }

  return false;
}

bool CollisionDetector::appendCollision(EdgeProxy* edge_a,
                                        const FaceProxy* triangle_b) {
  //    std::cout << "Broad phase hit: " << *edge_a << " vs " << *triangle_b <<
  //    '\n';

  std::list<ContinuousTimeCollision*> potentialCollisions;

  if (triangle_b->allApicesEnabled()) {
    VertexFaceCollision vf1(edge_a->getStrandPointer(),
                            edge_a->getVertexIndex() + 1, triangle_b);
    if (vf1.analyse())
      potentialCollisions.push_back(new VertexFaceCollision(vf1));

    VertexFaceCollision vf2(edge_a->getStrandPointer(),
                            edge_a->getVertexIndex(), triangle_b);
    if (vf2.analyse())
      potentialCollisions.push_back(new VertexFaceCollision(vf2));
  }

  for (int side = 0; side < 3; ++side) {
    const short side_p = (side + 1) % 3;

    if (triangle_b->hasEnabledApex(side) &&
        triangle_b->hasEnabledApex(side_p)) {
      EdgeFaceCollision ef(edge_a->getStrandPointer(), edge_a->getVertexIndex(),
                           triangle_b, triangle_b->getVertexIdx(side),
                           triangle_b->getVertexIdx(side_p), side, side_p);
      if (ef.analyse())
        potentialCollisions.push_back(new EdgeFaceCollision(ef));
    }
  }

  bool atLeastOnePositive = !potentialCollisions.empty();

  if (atLeastOnePositive) {
#pragma omp critical(pushCTCollision)
    {
      m_continuousTimeCollisions.insert(m_continuousTimeCollisions.end(),
                                        potentialCollisions.begin(),
                                        potentialCollisions.end());
    }
  }

  return atLeastOnePositive;
}

template <>
bool CollisionDetector::appendCollision<EdgeFaceIntersection>(
    ElementProxy* elem_a, ElementProxy* elem_b) {
  EdgeProxy* const edge_a = dynamic_cast<EdgeProxy*>(elem_a);

  if (edge_a) {
    FaceProxy* const triangle_b = dynamic_cast<FaceProxy*>(elem_b);
    return triangle_b && appendIntersection(edge_a, triangle_b);
  } else {
    EdgeProxy* const edge_b = dynamic_cast<EdgeProxy*>(elem_b);
    if (edge_b) {
      FaceProxy* const triangle_a = dynamic_cast<FaceProxy*>(elem_a);
      return triangle_a && appendIntersection(edge_b, triangle_a);
    }
  }
  return false;
}

bool CollisionDetector::appendIntersection(EdgeProxy* edge_a,
                                           const FaceProxy* triangle_b) {
  if (!triangle_b->allApicesEnabled()) {
    return false;
  }

  EdgeFaceIntersection intersection(
      edge_a->getStrandPointer(), edge_a->getVertexIndex(), edge_a, triangle_b);

  if (intersection.analyse()) {
    EdgeFaceIntersection* pInt = new EdgeFaceIntersection(intersection);
#pragma omp critical(pushProxCollision)
    { m_proximityCollisions.push_back(pInt); }
    return true;
  }

  return false;
}

void CollisionDetector::computeCollisions(const BVHNodeType& node_a,
                                          const BVHNodeType& node_b) {
  // If the bounding volumes do not overlap, there are no possible collisions
  // between their objects
  if (!intersect(node_a.BBox(), node_b.BBox())) return;

  // If both bounding volumes are leaves, add their contents to list potential
  // collisions
  if (node_a.IsLeaf() && node_b.IsLeaf()) {
    if (&node_a ==
        &node_b)  // As DK noticed, this case still needs to be considered. We
                  // still can avoid colliding an element with itself though.
    {
      const uint32_t leaf_a_begin = node_a.LeafBegin();
      const uint32_t leaf_a_end = node_a.LeafEnd();

      if (!m_ignoreContinuousTime) {
        for (unsigned int i = leaf_a_begin; i < leaf_a_end; ++i)
          for (unsigned int j = leaf_a_begin; j < i; ++j) {
            appendCollision<ContinuousTimeCollision>(m_elementProxies[i],
                                                     m_elementProxies[j]);
          }
      }
      if (!m_ignoreProximity) {
        for (unsigned int i = leaf_a_begin; i < leaf_a_end; ++i)
          for (unsigned int j = leaf_a_begin; j < i; ++j) {
            appendCollision<EdgeFaceIntersection>(m_elementProxies[i],
                                                  m_elementProxies[j]);
          }
      }
    } else {
      const uint32_t leaf_a_begin = node_a.LeafBegin();
      const uint32_t leaf_a_end = node_a.LeafEnd();
      const uint32_t leaf_b_begin = node_b.LeafBegin();
      const uint32_t leaf_b_end = node_b.LeafEnd();

      if (!m_ignoreContinuousTime) {
        for (unsigned int i = leaf_a_begin; i < leaf_a_end; ++i)
          for (unsigned int j = leaf_b_begin; j < leaf_b_end; ++j) {
            appendCollision<ContinuousTimeCollision>(m_elementProxies[i],
                                                     m_elementProxies[j]);
          }
      }
      if (!m_ignoreProximity) {
        for (unsigned int i = leaf_a_begin; i < leaf_a_end; ++i)
          for (unsigned int j = leaf_b_begin; j < leaf_b_end; ++j) {
            appendCollision<EdgeFaceIntersection>(m_elementProxies[i],
                                                  m_elementProxies[j]);
          }
      }
    }
  }

  // If one bounding volume is a leaf, we must recurse on the other volume
  else if (node_a.IsLeaf()) {
    computeCollisions(node_a, m_bvh.GetNode(node_b.ChildIndex()));
    computeCollisions(node_a, m_bvh.GetNode(node_b.ChildIndex() + 1));
  } else if (node_b.IsLeaf()) {
    computeCollisions(m_bvh.GetNode(node_a.ChildIndex()), node_b);
    computeCollisions(m_bvh.GetNode(node_a.ChildIndex() + 1), node_b);
  } else {
    computeCollisions(m_bvh.GetNode(node_a.ChildIndex()),
                      m_bvh.GetNode(node_b.ChildIndex()));
    computeCollisions(m_bvh.GetNode(node_a.ChildIndex() + 1),
                      m_bvh.GetNode(node_b.ChildIndex()));
    if (&node_a != &node_b)  // We need only to explore one side of the diagonal
      computeCollisions(m_bvh.GetNode(node_a.ChildIndex()),
                        m_bvh.GetNode(node_b.ChildIndex() + 1));
    computeCollisions(m_bvh.GetNode(node_a.ChildIndex() + 1),
                      m_bvh.GetNode(node_b.ChildIndex() + 1));
  }
}

void CollisionDetector::filterWithSpatialHashMap(const Scalar largestBBox) {
  if (!m_sortedAABBFunctor) {
    m_sortedAABBFunctor = new ElementProxySortedAABBFunctor(m_elementProxies);
  }
  if (!m_hashMap) {
    m_hashMap = new SpatialHashMapT();
  }

  m_hashMap->setCellSize((5 * largestBBox) > 1.0 ? (5 * largestBBox) : 1.0);

  const ElementProxySortedAABBFunctor& func(*m_sortedAABBFunctor);

  const unsigned numEdgeProxies = func.numEdgeProxies();

  //    Timer tt( "CollisionDetector::filterWithSpatialHashMap" ) ;
  m_hashMap->addToFootPrint(func, 0, numEdgeProxies);

#pragma omp parallel for
  for (int k = numEdgeProxies; k < m_elementProxies.size(); ++k) {
    if (!m_hashMap->isInFootPrint(func, k, k + 1)) {
      func.elementProxy(k)->resetBoundingBox();
    }
  }

  //    m_hashMap->clear() ;
  delete m_hashMap;
  m_hashMap = NULL;
}

} /* namespace strandsim */
