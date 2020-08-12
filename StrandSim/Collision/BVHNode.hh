/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TREES_BVHNODE_HH_S
#define TREES_BVHNODE_HH_S

namespace strandsim {

template <typename BBoxT>
struct BVHNode {
  typedef BBoxT BBoxType;
  typedef typename BBoxT::PointType PointType;

  BVHNode() : m_end((unsigned int)(-1)) {}

  BVHNode(const BBoxType& bbox, const unsigned int child_index)
      : m_bbox(bbox), m_index(child_index), m_end((unsigned int)(-1)) {}

  BVHNode(const BBoxType& bbox, const unsigned int leaf_begin,
          const unsigned int leaf_end)
      : m_bbox(bbox), m_index(leaf_begin), m_end(leaf_end) {}

  bool IsLeaf() const { return m_end != (unsigned int)(-1); }

  unsigned int LeafBegin() const { return m_index; }

  unsigned int LeafEnd() const { return m_end; }

  unsigned int ChildIndex() const { return m_index; }

  const BBoxType& BBox() const { return m_bbox; }

  BBoxType& BBox() { return m_bbox; }

  void SetBBox(const BBoxType& bbox) { m_bbox = bbox; }

  BBoxType m_bbox;       ///< the node's bbox
  unsigned int m_index;  ///< if an INNER node: the node's first child index, if
                         ///< LEAF: the leaf begin index
  unsigned int m_end;  ///< if an INNER node: -1, if a LEAF: the leaf end index
};

}  // namespace strandsim

#endif
