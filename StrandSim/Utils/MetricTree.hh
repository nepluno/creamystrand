/**
 * \copyright 2012 Gilles Daviet
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef METRICTREE_HH
#define METRICTREE_HH

#include <algorithm>
#include <limits>
#include <list>
#include <queue>
#include <vector>

#ifdef SS_MTREE_VERBOSE
#define SS_MTREE_LOG_LINE(a)     \
  do {                           \
    std::cerr << a << std::endl; \
  } while (0)
#else
#define SS_MTREE_LOG_LINE(a)
#endif

namespace strandsim {

/*!
  \brief Acceleration structure for non-euclidian metrics.
  Works by constructing a hierarchical structure of bounding spheres
  and using the triangular inequality to prune out branches at query time.
  Based on a classical M-Tree with some ideas from Sparse Spatial Selection [
  Brisaboa et al 2008 ] \tparam ScalarT Type of the codomain of the metric (
  usually \c float or \c double ) \tparam DataT Type of the arguments ( data
  points ) of the metric. \tparam LeafDataT Tag that can be associated with the
  objects inserted in the tree. Those tags will be retrieved when running
  queries on the MetricTree. \taparam MetricT Type of the metric functor The
  following signature should be defined: ScalarT MetricT::operator() (const
  Data&, const DataT& ) ;
*/
template <typename ScalarT, typename DataT, typename LeafDataT,
          typename MetricT>
class MetricTree {
 public:
  //! Index of an entry ; useful for easily retrieving objects already inserted
  //! in the tree.
  typedef size_t EntryIndexT;

  //! Maximum number of children for one node of the tree. Thie higher M, the
  //! flatter the tree
  static const int M = 24;
  //! When trying to insert a object into a node N with existing children, if
  //! its distance to the center of all existing children is greater than
  //! alpha()*diam( N ), create a new node
  inline static ScalarT alpha() { return 0.275; }

  //! Constructor
  /*! \param metric The functor used as the metric
    \param nObj Approximate number of objects that will be inserted, for
    performance purposes \param Approximate diameter of the objects that will be
    inserted, for performance purposes
    */
  MetricTree(const MetricT& metric, size_t nObj, ScalarT diam)
      : m_metric(metric), m_diam(diam), m_nCalls(0) {
    m_leafData.reserve(nObj);
    m_leafNodeIndexes.reserve(nObj);
    m_entries.reserve(2 * nObj);
    m_nodes.reserve(m_entries.size() / M);

    clear();
  }

  ~MetricTree() {}

  //! Clears the tree
  void clear() {
    m_leafData.clear();
    m_entries.clear();
    m_nodes.clear();

    m_nCalls = 0;
    m_size = 0;
    m_rootIndex = newNode();
  }

  //! Number of data points inside the tree
  size_t size() const { return m_size; }

  //! Inserts an data \p point with the tag \p leafData
  /*! \returns the entry index of the inserted object */
  EntryIndexT insert(const DataT& point, const LeafDataT& leafData) {
    EntryIndexT eIdx = m_entries.size();
    m_entries.push_back(Entry(point, m_leafData.size()));
    m_leafData.push_back(leafData);

    m_leafNodeIndexes.push_back(insert(m_rootIndex, point, eIdx, m_diam));

    m_size++;

    return eIdx;
  }

  //! Updates the position of a data point
  void update(const EntryIndexT eIdx, const DataT& newPoint) {
    assert(eIdx < m_entries.size());

    Entry& entry = m_entries[eIdx];
    const NodeIndexT nodeIdx = m_leafNodeIndexes[entry.leafData()];

    remove(nodeIdx, eIdx);

    entry.m_center = newPoint;
    m_leafNodeIndexes[entry.leafData()] =
        insert(m_rootIndex, newPoint, eIdx, m_diam);
  }

  //! Range query ; returns all the objects within a given \p radius of \p point
  /*! \p point does not need to have been inserted in the tree.
    \param leadfData The tags of all the objects of the tree that are with \p
    radius of \p point
    */
  void range(const DataT& point, const ScalarT radius,
             std::vector<LeafDataT>& leafData) const {
    std::list<QueryResult> indexes;
    range(m_rootIndex, point, radius, 0, indexes);

    leafData.reserve(indexes.size());
    for (auto it = indexes.begin(); it != indexes.end(); ++it) {
      leafData.push_back(m_leafData[it->m_leafDataIndex]);
    }
  }

  //! kNearestNeighbours query; finds the \p k entries closest to the \p point
  //! and return their tags
  void kNN(const DataT& point, const unsigned k,
           std::vector<LeafDataT>& leafData) const {
    std::priority_queue<QueryResult> indexes;
    ScalarT radius = std::numeric_limits<ScalarT>::infinity();
    kNN(k, m_rootIndex, point, radius, 0, indexes);

    assert(indexes.size() == k);

    leafData.resize(indexes.size());
    for (int i = indexes.size(); i-- > 0;) {
      leafData[i] = m_leafData[indexes.top().m_leafDataIndex];
      indexes.pop();
    }
  }

  //! kNearestNeighbours query; finds the \p k entries closest to the \p point
  //! and return their tags
  /*! Implementation as successive range-queries with increasing radii ; for
   * comparsion purposes only
   */
  void old_kNN(const DataT& point, unsigned k,
               std::vector<LeafDataT>& leafData) const {
    Scalar radius = 1.e-2;

    k = std::min((unsigned)size(), k);

    std::list<QueryResult> indexes;
    do {
      indexes.clear();
      range(m_rootIndex, point, radius, 0, indexes);
      radius *= 10;
    } while (indexes.size() < k);

    indexes.sort();

    leafData.reserve(k);
    auto indexIt = indexes.begin();
    for (unsigned i = 0; i < k; ++i) {
      leafData.push_back(m_leafData[(indexIt++)->m_leafDataIndex]);
    }
  }

  //! Returns the tag associated with the data point that has entry index \p
  //! eIdx
  const LeafDataT& leafData(const EntryIndexT eIdx) {
    assert(eIdx < m_entries.size());
    assert(m_entries[eIdx].isLeaf());

    return m_leafData[m_entries[eIdx].leafData()];
  }

  //! Returns the data point with the entry index \p eIdx
  const DataT& point(const EntryIndexT eIdx) {
    assert(eIdx < m_entries.size());

    return m_entries[eIdx].m_center;
  }

  //! Displays a flat overview of the tree structure
  std::ostream& print(std::ostream& os) const {
    // TODO better formating
    os << " === Nodes === " << std::endl;
    for (unsigned i = 0; i < m_nodes.size(); ++i) {
      const Node& node = m_nodes[i];
      if (isRoot(node))
        os << " * ";
      else
        os << "   ";
      os << i << "\t";
      os << node.m_parent << "\t[";
      for (unsigned k = 0; k < node.m_size; ++k) {
        os << node.m_entries[k] << "\t";
      }
      os << "]" << std::endl;
    }

    os << " === Entries === " << std::endl;
    for (unsigned i = 0; i < m_entries.size(); ++i) {
      const Entry& entry = m_entries[i];
      if (entry.isLeaf())
        os << " * ";
      else
        os << "   ";
      os << i << "\t";
      if (entry.isLeaf())
        os << m_leafData[entry.leafData()];
      else
        os << entry.m_node;
      os << "\t";
      os << entry.m_dist << "\t";
      os << entry.m_radius << std::endl;
      os << "   \t" << entry.m_center << std::endl;
    }

    os << " ====== " << std::endl;
    os << " Population: " << m_size << " / Distance calls: " << m_nCalls
       << std::endl;

    return os;
  }

  //! Returns the total number of calls to the metric fuction ( for stats )
  void resetNCalls() { m_nCalls = 0; }

  //! Displays a few stats about the tree
  void printStats() const {
    std::cout << " Population: " << m_size << " / Distance calls: " << m_nCalls
              << std::endl;
    std::cout << " Num entries: " << m_entries.size()
              << " / Num nodes: " << m_nodes.size() << std::endl;
    std::cout << " Entries per node: "
              << (((ScalarT)m_entries.size()) / m_nodes.size()) << std::endl;
    std::cout << " Entries per obj: " << (((ScalarT)m_entries.size()) / m_size)
              << std::endl;
  }

 private:
  typedef int32_t NodeIndexT;
  typedef int32_t LeafDataIndexT;
  struct QueryResult {
    LeafDataIndexT m_leafDataIndex;
    ScalarT m_distance;

    QueryResult() {}

    QueryResult(LeafDataIndexT ldi, ScalarT d)
        : m_leafDataIndex(ldi), m_distance(d) {}

    bool operator<(const QueryResult& rhs) const {
      return m_distance < rhs.m_distance;
    }
  };
  struct PrioritizedEntry {
    EntryIndexT m_entryIndex;
    ScalarT m_distance;

    bool operator<(const PrioritizedEntry& rhs) const {
      return m_distance < rhs.m_distance;
    }
  };

  //! Tree node
  /*!
    A node is a list of pointer to children entries.
    Each non-leaf entry contains a pointer to a child node.
    */
  class Node {
   public:
    //! Number of children
    short m_size;
    //! Index of parent node
    NodeIndexT m_parent;
    //! Indexes of chil entries
    EntryIndexT m_entries[M];

    Node() : m_size(0), m_parent(-1) {}

    bool full() const { return m_size == M; }
    bool empty() const { return m_size == 0; }

    void add(EntryIndexT e) { m_entries[m_size++] = e; }

    void clear() { m_size = 0; }
  };

  //! Entry: Metric sphere ( or point if representing a leaf )
  class Entry {
   public:
    Entry() : m_dist(0.), m_radius(0.), m_node(-1) {}

    Entry(const DataT& point, LeafDataIndexT ldIdx)
        : m_dist(0), m_radius(0.), m_center(point), m_node(-1 - ldIdx) {}

    bool isLeaf() const { return m_node < 0; }
    NodeIndexT leafData() const { return -1 - m_node; }

    //! Distance to the center its parent entry
    ScalarT m_dist;
    //! Radius of the sphere
    ScalarT m_radius;
    //! Center of the sphere
    DataT m_center;
    //! Pointer to child node
    NodeIndexT m_node;
  };

  bool isLeaf(const Node& node) const {
    return node.empty() || m_entries[node.m_entries[0]].isLeaf();
  }

  bool isRoot(const Node& node) const { return node.m_parent == -1; }

  //! Wrapper above the metric in order to count calls
  ScalarT distance(const DataT& P1, const DataT& P2) const {
    //#pragma omp atomic // We do not really need an accurate count ( except for
    //profiling )
    ++m_nCalls;
    return m_metric(P1, P2);
  }

  //! Returns the parent entry of a given node
  /*! Need to scan all the children of the parent node */
  EntryIndexT parentEntry(NodeIndexT nodeIdx) const {
    assert(!isRoot(m_nodes[nodeIdx]));

    const Node& parent = m_nodes[m_nodes[nodeIdx].m_parent];

    // We should already know that, since we did insert top-down.
    // We could probably pass the parent index as a parameter, or
    // analyse the stack... that would be fun
    for (short k = 0; k < parent.m_size; ++k) {
      if (m_entries[parent.m_entries[k]].m_node == nodeIdx) {
        return parent.m_entries[k];
      }
    }

    //        std::cerr << " Did not find parent for node #" << nodeIdx <<
    //        std::endl ;
    assert(false);
    return -1;
  }

  //! Re-assign the children of entries \param pe1 and \param pe2,
  /*! eventually chosing new centers for those entries at the same time.
    \param createRoot If true, we have just created a new root and those entries
    are directly below it. Therefore, we cannot use the center of the parent
    entry to help us find a good partition
    */
  void partition(Entry& pe1, Entry& pe2, const bool createRoot) {
    // Dispatch child entries between e1 and e2

    Node& n1 = m_nodes[pe1.m_node];
    Node& n2 = m_nodes[pe2.m_node];

    // 0 - Construct full list of entries and clear nodes
    EntryIndexT list[M + 1];
    std::memcpy(list, n1.m_entries, M * sizeof(EntryIndexT));
    list[M] = n2.m_entries[0];

    n1.clear();
    n2.clear();

    ScalarT d1[M + 1];
    ScalarT d2[M + 1];

    // I - chose two entries to promote

    if (createRoot) {
      // a\ We've just created a new root -- we have no parents to compare to
      // Use the "least-converage" heuristic

      ScalarT ds[(M + 1) * (M + 1)];

      // Compute all distances ; we can afford to do that for the root
      for (short i = 0; i < M + 1; ++i) {
        ds[i + (M + 1) * i] = 0;
        for (short j = i + 1; j < M + 1; ++j) {
          ds[i + (M + 1) * j] = ds[j + (M + 1) * i] = distance(
              m_entries[list[i]].m_center, m_entries[list[j]].m_center);
        }
      }

      ScalarT minCov = std::numeric_limits<Scalar>::infinity();
      short mink1, mink2;

      for (short i = 0; i < M + 1; ++i) {
        for (short j = i + 1; j < M + 1; ++j) {
          ScalarT maxRad = 0;

          const ScalarT* di = ds + (M + 1) * i;
          const ScalarT* dj = ds + (M + 1) * j;

          for (unsigned k = 0; k < M + 1; ++k) {
            maxRad = std::max(maxRad, std::min(di[k], dj[k]));
          }

          //                    std::cout << " Max rad for # " << i << " " << j
          //                    << " -> " << maxRad << std::endl ;

          if (maxRad < minCov) {
            minCov = maxRad;
            mink1 = i;
            mink2 = j;
          }
        }
      }

      pe1.m_center = m_entries[list[mink1]].m_center;
      pe2.m_center = m_entries[list[mink2]].m_center;

      memcpy(d1, ds + (M + 1) * mink1, (M + 1) * sizeof(ScalarT));
      memcpy(d2, ds + (M + 1) * mink2, (M + 1) * sizeof(ScalarT));

      SS_MTREE_LOG_LINE(" Re-part root with children " << mink1 << " and "
                                                       << mink2);

    } else {
      //   b\ Parent already exist -- find the furthest entry

      ScalarT maxd = -1;
      short maxk = -1;

      for (short k = 0; k < M + 1; ++k) {
        d1[k] = m_entries[list[k]].m_dist;
        if (d1[k] > maxd) {
          maxd = d1[k];
          maxk = k;
        }
      }

      pe2.m_center = m_entries[list[maxk]].m_center;
      pe2.m_radius = m_entries[list[maxk]].m_radius;

      SS_MTREE_LOG_LINE(" Partionning against child " << maxk);

      for (short k = 0; k < M + 1; ++k) {
        if (k == maxk) {
          d2[k] = 0.;
        } else {
          Entry& e = m_entries[list[k]];
          d2[k] = distance(pe2.m_center, e.m_center);
        }
      }

      // Oops, we forgot to update the distance of new entries to their parent
      if (!isRoot(m_nodes[n1.m_parent])) {
        const EntryIndexT gpeIdx = parentEntry(n1.m_parent);
        const Entry& gpe = m_entries[gpeIdx];
        pe2.m_dist = distance(pe2.m_center, gpe.m_center);
      }
    }

    ScalarT oldRadius = pe1.m_radius + pe2.m_radius;

    pe1.m_radius = 0.;
    pe2.m_radius = 0.;

    // II - Partition -- Use generalized hyperplane policy

    for (short k = 0; k < M + 1; ++k) {
      Entry& e = m_entries[list[k]];
      Entry& pe = d1[k] < d2[k] ? pe1 : pe2;
      Node& n = m_nodes[pe.m_node];

      n.add(list[k]);
      e.m_dist = std::min(d1[k], d2[k]);  // FIXME
      pe.m_radius = std::max(pe.m_radius, e.m_dist + e.m_radius);

      if (!e.isLeaf()) {
        m_nodes[e.m_node].m_parent = pe.m_node;
      }
    }

    // Everything shoud be alright

    if (!createRoot) {
      pe1.m_radius = std::min(oldRadius, pe1.m_radius);
      //            pe2.m_radius = std::min( oldRadius, pe2.m_radius ) ;
      SS_MTREE_LOG_LINE("we should not increase this " << oldRadius);
    }
  }

  NodeIndexT newNode() {
    const NodeIndexT nodeIdx = m_nodes.size();
    m_nodes.push_back(Node());

    return nodeIdx;
  }

  EntryIndexT newEntry() {
    const EntryIndexT entryIdx = m_entries.size();
    m_entries.push_back(Entry());

    return entryIdx;
  }

  void createParentNodeAndEntry(NodeIndexT& newNodeIdx,
                                EntryIndexT& newEntryIdx,
                                const EntryIndexT childEntryIdx) {
    newNodeIdx = newNode();
    newEntryIdx = newEntry();

    m_nodes.back().add(childEntryIdx);
    m_entries.back().m_node = newNodeIdx;

    Entry& child = m_entries[childEntryIdx];
    if (!child.isLeaf()) {
      m_nodes[child.m_node].m_parent = newNodeIdx;
    }
  }

  //! Splits node \param nodeIdx to accomodate the entry \param eIdx
  NodeIndexT split(const NodeIndexT nodeIdx, const EntryIndexT eIdx) {
    SS_MTREE_LOG_LINE(" Splitting node #" << nodeIdx << " to add entry #"
                                          << eIdx);

    NodeIndexT parentIdx;
    EntryIndexT peIdx;

    bool createRoot;
    if ((createRoot = isRoot(m_nodes[nodeIdx])))  // Pretty sure this is meant
                                                  // as an assignment :)
    {
      m_rootIndex = newNode();

      SS_MTREE_LOG_LINE(" Creating new root at index " << m_rootIndex);

      parentIdx = m_rootIndex;
      m_nodes[nodeIdx].m_parent = parentIdx;

      peIdx = newEntry();

      m_entries.back().m_node = nodeIdx;
      m_nodes[m_rootIndex].add(peIdx);
    } else {
      parentIdx = m_nodes[nodeIdx].m_parent;
      peIdx = parentEntry(nodeIdx);
    }

    NodeIndexT newNodeIdx;
    EntryIndexT newEntryIdx;
    createParentNodeAndEntry(newNodeIdx, newEntryIdx, eIdx);

    Entry& e1 = m_entries[peIdx];
    Entry& e2 = m_entries[newEntryIdx];

    partition(e1, e2, createRoot);

    Node& parent = m_nodes[parentIdx];

    if (parent.full()) {
      return split(parentIdx, newEntryIdx);
    } else {
      parent.add(newEntryIdx);
      m_nodes[newNodeIdx].m_parent = parentIdx;

      return parentIdx;
    }
  }

  //! Insertes a new entry in the tree
  /*!
    \param nodeIdx The node in which to insert the entry
    \param eIdx Pointer to then entry to insert
    \param point Center of the entry to insert
    \param diam Diameter of nodeIdx
    \param dist Distance of \p point to the center of the parent entry of \p
    nodeIdx
    */
  NodeIndexT insert(const NodeIndexT nodeIdx, const DataT& point,
                    const EntryIndexT eIdx, const ScalarT diam,
                    const ScalarT dist = 0.) {
    SS_MTREE_LOG_LINE(" Insert " << eIdx << " in node #" << nodeIdx
                                 << ", dist = " << dist);

    Node& node = m_nodes[nodeIdx];

    if (isLeaf(node)) {
      m_entries[eIdx].m_dist = dist;

      if (node.full()) {
        SS_MTREE_LOG_LINE(" Node is full leaf");

        return split(nodeIdx, eIdx);

      } else {
        node.add(eIdx);

        SS_MTREE_LOG_LINE(" Node is non-full leaf, inserting entry #" << eIdx);
      }
      return nodeIdx;

    } else {
      ScalarT d[M];
      short mink = -1;
      ScalarT mind = std::numeric_limits<ScalarT>::infinity();

      // TODO We could probably use the dist parameter  + triangular inequality
      // to avoid some distance compuations

      for (short k = 0; k < node.m_size; ++k) {
        const Entry& e = m_entries[node.m_entries[k]];

        const Scalar lmd = std::min(e.m_radius, mind);

        if (lmd + dist < e.m_dist) {
          d[k] = 0;
        } else {
          d[k] = distance(point, e.m_center);
          // Closest satisfyning center
          if (d[k] <= lmd) {
            mink = k;
            mind = d[k];
          }
        }
      }

      if (mink == -1) {
        ScalarT minInc = std::numeric_limits<ScalarT>::infinity();
        for (short k = 0; k < node.m_size; ++k) {
          const Entry& e = m_entries[node.m_entries[k]];

          if (minInc + dist < e.m_dist) continue;

          if (d[k] == 0) d[k] = distance(point, e.m_center);
          const ScalarT inc = d[k];  // - e.m_radius ;

          if (inc < minInc) {
            minInc = inc;
            mink = k;
          }
        }

        assert(mink != -1);

        if (d[mink] > alpha() * diam) {
          SS_MTREE_LOG_LINE(" Too far away, adding node ( " << d[mink] << " vs "
                                                            << diam);

          NodeIndexT newNodeIdx;
          EntryIndexT newEntryIdx;
          createParentNodeAndEntry(newNodeIdx, newEntryIdx, eIdx);

          m_entries[eIdx].m_dist = 0.;
          m_entries[newEntryIdx].m_dist = dist;
          m_entries[newEntryIdx].m_center = point;

          m_nodes[newNodeIdx].m_parent = nodeIdx;

          Node& node = m_nodes[nodeIdx];
          if (node.full()) {
            return split(nodeIdx, newEntryIdx);
          } else {
            node.add(newEntryIdx);
            return newNodeIdx;
          }

        } else {
          SS_MTREE_LOG_LINE(" No satisfying entry, expanding child " << mink);

          m_entries[node.m_entries[mink]].m_radius = d[mink];

          return insert(m_entries[node.m_entries[mink]].m_node, point, eIdx,
                        2. * d[mink], d[mink]);
        }
      } else {
        SS_MTREE_LOG_LINE(" Inserting in satisfying child " << mink);

        return insert(m_entries[node.m_entries[mink]].m_node, point, eIdx,
                      2. * m_entries[node.m_entries[mink]].m_radius, d[mink]);
      }
    }
  }

  //! Range query
  /*!
    \param nodeIdx Current node to examine
    \param point Center of the query
    \param radius Radius of the query
    \param dist Distance of \p point to the center of the parent entry of \p
    nodeIdx
    */
  void range(const NodeIndexT nodeIdx, const DataT& point, const ScalarT radius,
             const ScalarT dist, std::list<QueryResult>& indexes) const {
    const Node& node = m_nodes[nodeIdx];

    for (short k = 0; k < node.m_size; ++k) {
      const Entry& e = m_entries[node.m_entries[k]];

      if (e.m_radius + radius < std::fabs(dist - e.m_dist)) {
        continue;
      }

      const ScalarT d = distance(point, e.m_center);

      if (d <= radius + e.m_radius) {
        if (e.isLeaf()) {
          indexes.push_back(QueryResult(e.leafData(), d));
        } else {
          range(e.m_node, point, radius, d, indexes);
        }
      }
    }
  }

  //! k-Nearest-Neighbours query
  /*!
    \param K the number or nearest-neighbours to retrieve
    \param nodeIdx Current node to examine
    \param point Center of the query
    \param radius Current radius of the query ; will be refined as more nearest
    neighbours are found \param dist Distance of \p point to the center of the
    parent entry of \p nodeIdx
    */
  void kNN(const unsigned K, const NodeIndexT nodeIdx, const DataT& point,
           ScalarT& radius, const ScalarT dist,
           std::priority_queue<QueryResult>& indexes) const {
    const Node& node = m_nodes[nodeIdx];

    PrioritizedEntry sortedEntries[M];
    short nEntries = 0;

    // Sort entries by distance to query point
    for (short k = 0; k < node.m_size; ++k) {
      const Entry& e = m_entries[node.m_entries[k]];

      if (e.m_radius + radius < std::fabs(dist - e.m_dist)) {
        continue;
      }

      const ScalarT d = distance(point, e.m_center);

      const PrioritizedEntry priorityEntry = {node.m_entries[k], d};
      sortedEntries[nEntries++] = priorityEntry;
    }

    std::sort(sortedEntries, sortedEntries + nEntries);

    // Depth-first recursive search on children
    for (short k = 0; k < nEntries; ++k) {
      const PrioritizedEntry& priorityEntry = sortedEntries[k];

      const Entry& e = m_entries[priorityEntry.m_entryIndex];
      const ScalarT d = priorityEntry.m_distance;

      // We recheck this inequality as the query radius may have become smaller
      if (e.m_radius + radius < std::fabs(dist - e.m_dist)) {
        continue;
      }

      if (d <= radius + e.m_radius) {
        if (e.isLeaf()) {
          indexes.push(QueryResult(e.leafData(), d));
          if (indexes.size() > K) {
            indexes.pop();
            // Query radius refinement
            radius = indexes.top().m_distance;
          }
        } else {
          kNN(K, e.m_node, point, radius, d, indexes);
        }
      }
    }
  }

  //! Removes the entry \p eIdx fronm the node \p nodeIdx
  void remove(const NodeIndexT nodeIdx, const EntryIndexT eIdx) {
    Node& node = m_nodes[nodeIdx];

    int ke = -1;
    for (short k = 0; k < node.m_size; ++k) {
      if (node.m_entries[k] == eIdx) {
        ke = k;
        break;
      }
    }

    if (ke == -1) return;

    --node.m_size;
    std::swap(node.m_entries[ke], node.m_entries[node.m_size]);

    // TODO recomp radiuses ?
  }

  //! Number of data points in the tree
  size_t m_size;

  //! Index of the root of the tree in m_nodes
  NodeIndexT m_rootIndex;

  std::vector<Node> m_nodes;
  std::vector<Entry> m_entries;

  //! Tags data
  std::vector<LeafDataT> m_leafData;
  //! For each leaf, index of the node that immediately contains it
  std::vector<NodeIndexT> m_leafNodeIndexes;

  const MetricT& m_metric;
  const ScalarT m_diam;

  mutable unsigned m_nCalls;
};

template <typename ScalarT, typename DataT, typename LeafDataT,
          typename MetricT>
std::ostream& operator<<(
    std::ostream& os, const MetricTree<ScalarT, DataT, LeafDataT, MetricT>& M) {
  return M.print(os);
}

}  // namespace strandsim

#endif  // METRICTREE_HH
