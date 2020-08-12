/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef INDEXED_SET_HH_
#define INDEXED_SET_HH_

#include <map>
#include <set>
#include <vector>

template <typename IndexT, typename ElemT,
          typename ElemComparatorT = std::less<ElemT>>
struct IndexedSet {
  typedef std::set<ElemT, ElemComparatorT> SetType;

  enum { invalidIndex = -1 };

  IndexedSet() : m_index(invalidIndex) {}

  explicit IndexedSet(IndexT index) : m_index(index), m_set(SetType()) {}

  IndexT m_index;
  SetType m_set;
};

// A functor that will tell us if the IndexedSet's index already exists
template <typename IndexT>
struct IndexComparator {
  explicit IndexComparator(IndexT index) : m_index(index) {}

  template <typename ElemT>
  bool operator()(const IndexedSet<IndexT, ElemT>& swv) {
    return swv.m_index == m_index;
  }

 private:
  const IndexT m_index;
};

template <typename IndexT, typename ElemT,
          typename ElemComparatorT = std::less<ElemT>>
struct IndexedSetContainer {
  typedef IndexedSet<IndexT, ElemT, ElemComparatorT> IndexedSetType;
  typedef std::map<IndexT, IndexedSetType> ContainerImplType;

  typedef enum { nextSet } nextSetOnly;

  explicit IndexedSetContainer(
      const ContainerImplType& data = ContainerImplType())
      : m_impl(data) {}

  typedef struct const_iterator {
    typedef IndexedSet<IndexT, ElemT> IndexedSetType;

    const_iterator() {}

    const_iterator(typename ContainerImplType::const_iterator main,
                   typename IndexedSetType::SetType::const_iterator sub,
                   typename ContainerImplType::const_iterator mainEnd)
        : m_mainIndex(main), m_subIndex(sub), m_mainIndexEnd(mainEnd) {}

    // Constructpr initliazing subIndex to begin() iff mainIndex is valid
    const_iterator(typename ContainerImplType::const_iterator main,
                   typename ContainerImplType::const_iterator mainEnd)
        : m_mainIndex(main), m_mainIndexEnd(mainEnd) {
      if (main != mainEnd) {
        m_subIndex = m_mainIndex->second.m_set.begin();
      }
    }

    typename ContainerImplType::const_iterator m_mainIndex;
    typename IndexedSetType::SetType::const_iterator m_subIndex;

    typename ContainerImplType::const_iterator m_mainIndexEnd;

    const_iterator& operator++()  // prefix operator
    {
      m_subIndex++;
      if (m_subIndex == m_mainIndex->second.m_set.end()) {
        m_mainIndex++;
        if (m_mainIndex != m_mainIndexEnd) {
          m_subIndex = m_mainIndex->second.m_set.begin();
        }
      }
      return *this;
    }

    // This operator takes only one possible argument value (nextSet) and is
    // actually an increment operator for the strand.
    const_iterator& operator+=(int e) {
      m_mainIndex++;
      m_subIndex = m_mainIndex->second.m_set.begin();
      return *this;
    }

    std::pair<IndexT, ElemT> operator*() const {
      return std::pair<IndexT, ElemT>(m_mainIndex->first, *m_subIndex);
    }

    bool operator==(const const_iterator& other) {
      return (
          m_mainIndex == other.m_mainIndex &&
          (m_mainIndex == m_mainIndexEnd || m_subIndex == other.m_subIndex));
    }

    bool operator!=(const const_iterator& other) { return !(*this == other); }

    const typename IndexedSetType::SetType& getSet() const {
      return m_mainIndex->second.m_set;
    }

  } const_iterator;

  const_iterator begin() const {
    // This is a
    return const_iterator(m_impl.begin(), m_impl.end());
  }

  const_iterator end() const {
    return const_iterator(m_impl.end(), m_impl.end());
  }

  const_iterator find(const IndexT& index) const {
    return const_iterator(m_impl.find(index), m_impl.end());
  }

  const typename IndexedSetType::SetType& getSet(
      const const_iterator& swv) const {
    return swv.m_mainIndex->second.m_set;
  }

  void setSet(IndexT index, const typename IndexedSetType::SetType& set) {
    auto it = m_impl.find(index);
    if (it == m_impl.end()) {
      IndexedSetType& indexedSet = m_impl[index];
      indexedSet.m_index = index;
      indexedSet.m_set = set;
    } else {
      it->second.m_set = set;
    }
  }

  void insert(const std::pair<IndexT, ElemT>& toInsert) {
    auto it = m_impl.find(toInsert.first);
    if (it == m_impl.end()) {
      IndexedSetType& set = m_impl[toInsert.first];
      set.m_index = toInsert.first;
      set.m_set.insert(toInsert.second);
    } else {
      it->second.m_set.insert(toInsert.second);
    }
  }

  template <typename IteratorT>
  void insert(const IteratorT& begin, const IteratorT& end) {
    for (IteratorT it = begin; it != end; ++it) {
      insert(*it);
    }
  }

  template <typename IteratorT>
  void insert(IndexT index, const IteratorT& begin, const IteratorT& end) {
    auto it = m_impl.find(index);
    if (it == m_impl.end()) {
      IndexedSetType& indexedSet = m_impl[index];
      indexedSet.m_index = index;
      indexedSet.m_set.insert(begin, end);
    } else {
      it->second.m_set.insert(begin, end);
    }
  }

  void remove(const std::pair<IndexT, ElemT>& toRemove) {
    auto it = m_impl.find(toRemove.first);
    if (it != m_impl.end()) {
      typename IndexedSetType::SetType& foundSetImpl = it->second.m_set;
      foundSetImpl.erase(toRemove.second);
      if (foundSetImpl.empty()) m_impl.erase(it);
    }
  }

  // This function should not work on const_iterator
  // TODO implement mutable iterator class
  const_iterator erase(const const_iterator& it) {
    const_iterator next = it;
    ++next;
    remove(std::make_pair((*it).first, (*it).second));
    return next;
  }

  void toggle(const std::pair<IndexT, ElemT>& toToggle) {
    auto it = m_impl.find(toToggle.first);
    if (it == m_impl.end()) {
      IndexedSetType& set = m_impl[toToggle.first];
      set.m_index = toToggle.first;
      set.m_set.insert(toToggle.second);
    } else {
      typename IndexedSetType::SetType& foundSetImpl = it->second.m_set;
      foundSetImpl.erase(toToggle.second);
      if (foundSetImpl.empty()) m_impl.erase(it);
    }
  }

  void removeSet(const IndexT& toRemove) { m_impl.erase(toRemove); }

  void clear() { m_impl.clear(); }

  bool empty() const { return m_impl.empty(); }

  std::pair<IndexT, ElemT> front() const {
    assert(!empty());
    auto firstSet = begin();
    return std::make_pair((*firstSet).first, *firstSet.getSet().begin());
  }

  bool hasIndex(int index) { return m_impl.end() != m_impl.find(index); }

  // size is the number of sets, not the cardinal of their union
  typename ContainerImplType::size_type size() const { return m_impl.size(); }

  typename IndexedSetType::SetType::size_type unionCardinal() const {
    typename IndexedSetType::SetType::size_type sum = 0;
    auto add = [&sum](const std::pair<const int, IndexedSetType>& it) {
      sum += it.second.m_set.size();
    };
    std::for_each(m_impl.begin(), m_impl.end(), add);

    return sum;
  }

  template <typename VectorT>
  void getSetsIndices(VectorT& indices) const {
    indices.clear();
    indices.reserve(this->size());
    for (auto it = m_impl.begin(); it != m_impl.end(); ++it) {
      indices.push_back(it->first);
    }
  }

 private:
  ContainerImplType m_impl;
};

namespace strandsim {
typedef IndexedSetContainer<int, int> VerticesSet;
}

#endif /* INDEXED_SET_HH_ */
