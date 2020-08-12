/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "BVH.hh"

#include "ElementProxy.hh"

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
#endif
#include <boost/thread.hpp>

namespace strandsim {

template <typename BBoxFunctorT>
void BVHBuilder<BBoxFunctorT>::build(BBoxFunctorT& bboxes, BVH* bvh) {
  const uint32_t n = bboxes.size();

  const BBoxType root_bbox = parallel_compute_bbox(bboxes, 0u, n);
  m_stack.push_back(StackNode(0u,           // node index
                              0u,           // begin index
                              n,            // end index
                              0u,           // node depth
                              root_bbox));  // node bbox

  std::vector<BVHNodeType>& nodes = bvh->GetNodeVector();
  nodes.resize(1u);

  // Conservative estimate of the total number of nodes -- essential to avoid
  // data-races in parallelisation
  nodes.reserve(2 * n + 1);

  const unsigned max_threads = omp_get_max_threads();

  // First compute a root node per thread
  while (!m_stack.empty() && m_stack.size() < max_threads) {
    // FIFO ( to get good parallelism )
    StackNode node = m_stack.front();
    m_stack.pop_front();

    const BBoxType& node_bbox =
        node.m_begin + n == node.m_end
            ? root_bbox
            : parallel_compute_bbox(bboxes, node.m_begin, node.m_end);
    process_node(bboxes, node, node_bbox, m_stack, nodes);
  }

// Then parallel-compute the bvh from each of those root nodes
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < m_stack.size(); ++i) {
    StackNode root = m_stack[i];

    std::vector<BVHNodeType> local_nodes;
    local_nodes.reserve(2 * (root.m_end - root.m_begin) + 1);
    local_nodes.resize(1u);

    const uint32_t root_index = root.m_node_index;
    root.m_node_index = 0;

    std::deque<StackNode> stack;
    stack.push_back(root);

    while (!stack.empty()) {
      // LIFO ( to get branch locality )
      StackNode node = stack.back();
      stack.pop_back();

      const BBoxType& node_bbox =
          compute_bbox(bboxes, node.m_begin, node.m_end);
      process_node(bboxes, node, node_bbox, stack, local_nodes);
    }

    nodes[root_index] = local_nodes[0];

    if (!nodes[root_index].IsLeaf()) {
      // Adjust indices of thread-local nodes and insert them in global list
      unsigned offset;
#pragma omp critical
      {
        offset = nodes.size() - 1;
        nodes.resize(offset + local_nodes.size());
      }

      nodes[root_index].m_index += offset;

      // Adjust offset of inner nodes
      for (unsigned nidx = 1; nidx < local_nodes.size(); ++nidx) {
        if (!local_nodes[nidx].IsLeaf()) {
          local_nodes[nidx].m_index += offset;
        }
      }
      memcpy(&nodes[offset + 1], &local_nodes[1],
             (local_nodes.size() - 1) * sizeof(BVHNodeType));
    }
  }
}

template <typename BBoxFunctorT>
void BVHBuilder<BBoxFunctorT>::process_node(BBoxFunctorT& bboxes,
                                            const StackNode& node,
                                            const BBoxType& node_bbox,
                                            std::deque<StackNode>& stack,
                                            std::vector<BVHNodeType>& nodes) {
  if (node.m_end - node.m_begin <= m_max_leaf_size) {
    nodes[node.m_node_index] = BVHNodeType(node_bbox, node.m_begin, node.m_end);
    return;
  }

  const BBoxType& kd_bbox = presplit(node_bbox, node.m_kd_bbox);

  PointType edge = kd_bbox.max - kd_bbox.min;

  uint32_t split_dim =
      uint32_t(std::max_element(&edge[0], &edge[0] + 3) - &edge[0]);

  Scalar split_plane = (kd_bbox.max[split_dim] + kd_bbox.min[split_dim]) * 0.5;

  uint32_t split_index =
      partition(bboxes, node.m_begin, node.m_end, split_dim, split_plane);

  if (split_index == node.m_begin || split_index == node.m_end) {
    edge = node_bbox.max - node_bbox.min;
    split_dim = uint32_t(std::max_element(&edge[0], &edge[0] + 3) - &edge[0]);
    split_plane = (node_bbox.max[split_dim] + node_bbox.min[split_dim]) * 0.5;
    split_index =
        partition(bboxes, node.m_begin, node.m_end, split_dim, split_plane);
  }

  if (split_index == node.m_begin || split_index == node.m_end) {
    edge = node_bbox.max - node_bbox.min;
    split_dim = uint32_t(std::max_element(&edge[0], &edge[0] + 3) - &edge[0]);

    Scalar mean = 0.0;
    for (uint32_t i = node.m_begin; i < node.m_end; i++)
      mean += (bboxes[i].min[split_dim] + bboxes[i].max[split_dim]) * 0.5;

    split_plane = mean / Scalar(node.m_begin - node.m_end);
    split_index =
        partition(bboxes, node.m_begin, node.m_end, split_dim, split_plane);
  }

  if (split_index == node.m_begin || split_index == node.m_end) {
    nodes[node.m_node_index] = BVHNodeType(node_bbox, node.m_begin, node.m_end);
    return;
  }

  const uint32_t child_index = uint32_t(nodes.size());
  nodes.resize(child_index + 2u);

  BBoxType right_bbox = node.m_kd_bbox;
  right_bbox.min[split_dim] = split_plane;
  StackNode right_node(child_index + 1u, split_index, node.m_end,
                       node.m_depth + 1u, right_bbox);

  BBoxType left_bbox = node.m_kd_bbox;
  left_bbox.max[split_dim] = split_plane;
  StackNode left_node(child_index, node.m_begin, split_index, node.m_depth + 1u,
                      left_bbox);

  nodes[node.m_node_index] = BVHNodeType(node_bbox, child_index);

  stack.push_back(right_node);
  stack.push_back(left_node);
}

template <typename BBoxFunctorT>
BBoxType BVHBuilder<BBoxFunctorT>::presplit(const BBoxType& node_bbox,
                                            const BBoxType& kd_bbox) {
  int tests[3];
  for (uint32_t i = 0; i < 3; i++) tests[i] = 8;

  BBoxType out_bbox = kd_bbox;
  while (tests[0] && tests[1] && tests[2]) {
    const PointType& edge = out_bbox.max - out_bbox.min;
    const uint32_t split_dim =
        uint32_t(std::max_element(&edge[0], &edge[0] + 3) - &edge[0]);
    const Scalar split_plane =
        (out_bbox.min[split_dim] + out_bbox.max[split_dim]) * 0.5;

    tests[split_dim]--;
    if (split_plane < node_bbox.min[split_dim]) {
      out_bbox.min[split_dim] = split_plane;
      continue;
    }
    if (split_plane > node_bbox.max[split_dim]) {
      out_bbox.max[split_dim] = split_plane;
      continue;
    }
    return out_bbox;
  }
  return out_bbox;
}

void swap(BVH& a, BVH& b) { std::swap(a.GetNodeVector(), b.GetNodeVector()); }

template <typename BBoxFunctorT>
uint32_t partition(BBoxFunctorT& bboxes, const uint32_t begin,
                   const uint32_t end, const uint32_t axis, const float pivot) {
  uint32_t i = begin;
  uint32_t j = end - 1;
  while (i != j) {
    if (is_left(bboxes[i], axis, pivot) == false) {
      bboxes.swap(i, j);
      --j;
    } else
      ++i;
  }
  if (is_left(bboxes[i], axis, pivot) == true)
    return ++i;
  else
    return i;
}

template <typename BBoxFunctorT>
BBoxType parallel_compute_bbox(BBoxFunctorT& bboxes, const uint32_t begin,
                               const uint32_t end) {
  BBoxType node_bbox;

#pragma omp parallel
  {
    BBoxType loc_bbox;
#pragma omp for
    for (int i = begin; i < end; i++) loc_bbox.insert(bboxes[i]);
#pragma omp critical
    { node_bbox.insert(loc_bbox); }
  }

  return node_bbox;
}

template <typename BBoxFunctorT>
BBoxType compute_bbox(BBoxFunctorT& bboxes, const uint32_t begin,
                      const uint32_t end) {
  //    assert( begin != end ) ;

  BBoxType node_bbox(bboxes[begin]);

  for (uint32_t i = begin + 1; i < end; i++) node_bbox.insert(bboxes[i]);

  return node_bbox;
}

// Explicit template instantiations
template class BVHBuilder<ElementProxyBBoxFunctor>;

}  // namespace strandsim
