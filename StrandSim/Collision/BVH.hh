/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BVH_HH_S
#define BVH_HH_S

#include <algorithm>
#include <limits>
#include <stack>
#include <vector>

#include "../Core/BoundingBox.hh"
#include "BVHNode.hh"

namespace strandsim {

typedef BoundingBox<float> BBoxType;
typedef BVHNode<BBoxType> BVHNodeType;

class BVH {
 public:
  typedef BVHNodeType Node_Type;

  /// empty constructor
  BVH() {}

  /// returns the size of this object in bytes
  size_t ByteSize() const {
    return sizeof(BVHNodeType) * m_nodes.size() + sizeof(BVH);
  }

  /// get node vector
  const std::vector<BVHNodeType>& GetNodeVector() const { return m_nodes; }

  /// get node vector
  std::vector<BVHNodeType>& GetNodeVector() { return m_nodes; }

  /// get nodes pointer
  const BVHNodeType* GetNodes() const { return &m_nodes[0]; }

  /// get nodes pointer
  BVHNodeType* GetNodes() { return &m_nodes[0]; }

  /// get the i-th node
  const BVHNodeType& GetNode(const unsigned int i) const { return m_nodes[i]; }

  /// get the i-th node
  BVHNodeType& GetNode(const unsigned int i) { return m_nodes[i]; }

 private:
  std::vector<BVHNodeType> m_nodes;  ///< bvh nodes
};

void swap(BVH& a, BVH& b);

template <typename BBoxFunctorT>
class BVHBuilder {
 public:
  typedef BBoxType::PointType PointType;

  /// empty constructor
  BVHBuilder() : m_max_leaf_size(1u) {}

  void build(BBoxFunctorT& bboxes, BVH* bvh);

 private:
  BBoxType presplit(const BBoxType& node_bbox, const BBoxType& kd_bbox);

  struct StackNode {
    StackNode() {}
    StackNode(const unsigned int node, const unsigned int begin,
              const unsigned int end, const unsigned int depth,
              const BBoxType& kd_bbox)
        : m_node_index(node),
          m_begin(begin),
          m_end(end),
          m_depth(depth),
          m_kd_bbox(kd_bbox) {}

    unsigned int m_node_index;
    unsigned int m_begin;
    unsigned int m_end;
    unsigned int m_depth;
    BBoxType m_kd_bbox;
  };

  void process_node(BBoxFunctorT& bboxes, const StackNode& node,
                    const BBoxType& node_bbox, std::deque<StackNode>& stack,
                    std::vector<BVHNodeType>& nodes);

  BVH* m_bvh;                          ///< output bvh
  std::deque<StackNode> m_stack;       ///< internal stack
  const unsigned int m_max_leaf_size;  ///< maximum leaf size
};

template <typename BBoxFunctorT>
unsigned int partition(BBoxFunctorT& bboxes, const unsigned int begin,
                       const unsigned int end, const unsigned int axis,
                       const float pivot);

template <typename BBoxFunctorT>
BBoxType compute_bbox(BBoxFunctorT& bboxes, const unsigned int begin,
                      const unsigned int end);

template <typename BBoxFunctorT>
BBoxType parallel_compute_bbox(BBoxFunctorT& bboxes, const unsigned int begin,
                               const unsigned int end);
}  // namespace strandsim

#endif
