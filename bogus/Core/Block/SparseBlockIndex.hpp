/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_SPARSEBLOCKINDEX_HPP
#define BOGUS_SPARSEBLOCKINDEX_HPP

#include <algorithm>
#include <cassert>
#include <vector>

namespace bogus {

template <typename Derived>
struct SparseBlockIndexTraits {};

template <typename Derived>
struct SparseBlockIndexBase {
  typedef SparseBlockIndexTraits<Derived> Traits;
  typedef typename Traits::Index Index;
  typedef typename Traits::InnerIterator InnerIterator;

  //! Type of the array encoding the size of each block of the inner dimension.
  /*! which can be retrieved as \code innerOffsets[inner+1] -
   * innerOffsets[inner] \endcode */
  typedef std::vector<Index> InnerOffsetsType;

  //! Whether this index is currently valid
  bool valid;

  SparseBlockIndexBase(bool _valid = true) : valid(_valid) {}

  Derived& derived();
  const Derived& derived() const;

  //! Number of elemnts of the major indexing direction
  /*! i.e. number of rows for a row-major index */
  Index outerSize() const;
  //! Number of elements of the minor indexing direction
  /*! i.e. number of cols for a row-major index */
  Index innerSize() const;

  //! Returns whether the innerOffsetsArray() has been filled
  bool hasInnerOffsets() const;
  //! \sa InnerOFfsetsType
  const InnerOffsetsType& innerOffsetsArray() const;

  //! Returns the total number of nonZeros in this index
  /*! Depending on the index type, may performe som computations */
  Index nonZeros() const { return derived().nonZeros(); }

  //! Same as innerOffsetsArray, but returns a pointer instead. Assumes
  //! hasInnerOffsets()
  const Index* innerOffsetsData() const {
    assert(hasInnerOffsets());
    // innerOffsetsArray() guaranteed to not be empty, we can dereference its
    // first elt
    return &innerOffsetsArray()[0];
  }

  //! Returns the number of blocks at the outer index \p outerIdx
  Index size(const Index outerIdx) const;
  //! Returns an iterator to the first non-empty block of \p outerIdx
  InnerIterator begin(const Index outerIdx) const;
  //! Returns an iterator to the last non-empty block of \p outerIdx
  InnerIterator last(const Index outerIdx) const;
  //! Returns an iterator to the end of \p outerIdx
  InnerIterator end(const Index outerIdx) const;
};

//! Uncompressed sparse block index
template <bool Compressed, typename Index_, typename BlockPtr_ = Index_,
          template <typename> class ArrayType = ResizableSequenceContainer>
struct SparseBlockIndex
    : public SparseBlockIndexBase<
          SparseBlockIndex<Compressed, Index_, BlockPtr_, ArrayType> > {
  typedef Index_ Index;
  typedef BlockPtr_ BlockPtr;

  typedef SparseBlockIndexBase<
      SparseBlockIndex<Compressed, Index_, BlockPtr_, ArrayType> >
      Base;
  typedef typename Base::InnerOffsetsType InnerOffsetsType;
  typedef typename Base::InnerIterator InnerIterator;
  using Base::valid;

  //! Vector of ( inner index ; block pointer ) tuples encoding an inner vector
  typedef std::vector<std::pair<Index, BlockPtr> > Inner;
  //! Vector of inner vectors
  typedef typename ArrayType<Inner>::Type Outer;

  InnerOffsetsType innerOffsets;
  Outer outer;

  bool ordered;

  SparseBlockIndex() : Base(), ordered(true) {}

  void resizeOuter(Index size) { outer.resize(size); }
  void reserve(Index /*nnz*/) {}

  Index outerSize() const { return outer.size(); }
  const InnerOffsetsType& innerOffsetsArray() const { return innerOffsets; }

  template <bool Ordered>
  void insert(Index outIdx, Index inIdx, BlockPtr ptr) {
    outer[outIdx].push_back(std::make_pair(inIdx, ptr));
    ordered &= Ordered;
  }
  void insertBack(Index outIdx, Index inIdx, BlockPtr ptr) {
    insert<true>(outIdx, inIdx, ptr);
  }

  void finalize() {
    if (ordered) return;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
    for (int i = 0; i < (Index)outer.size(); ++i) {
      std::sort(outer[i].begin(), outer[i].end());
    }
  }

  void clear() {
    std::vector<Inner>(outer.size()).swap(outer);
    valid = true;
    ordered = true;
  }

  SparseBlockIndex& operator=(const SparseBlockIndex& o) {
    if (&o != this) {
      outer = o.outer;
      valid = o.valid;
      if (!o.innerOffsets.empty()) innerOffsets = o.innerOffsets;
    }
    return *this;
  }

  template <typename SourceDerived>
  SparseBlockIndex& operator=(
      const SparseBlockIndexBase<SourceDerived>& source);

  SparseBlockIndex& move(SparseBlockIndex& uncompressed) {
    if (&uncompressed != this) {
      outer.swap(uncompressed.outer);
      if (!uncompressed.innerOffsets.empty())
        innerOffsets.swap(uncompressed.innerOffsets);
      valid = uncompressed.valid;
      uncompressed.valid = false;
    }
    return *this;
  }

  template <typename SourceDerived>
  SparseBlockIndex& move(const SparseBlockIndexBase<SourceDerived>& source) {
    return (*this = source);
  }

  template <bool Symmetric, typename SourceDerived>
  SparseBlockIndex& setToTranspose(
      const SparseBlockIndexBase<SourceDerived>& source) {
    clear();
    resizeOuter(source.innerSize());
    valid = source.valid;

    for (typename SourceDerived::Index i = 0; i < source.outerSize(); ++i) {
      // For a symmetric matrix, do not store diagonal block in col-major index
      for (typename SourceDerived::InnerIterator it(source.derived(), i);
           it && (!Symmetric || it.inner() < i); ++it) {
        insertBack(it.inner(), i, it.ptr());
      }
    }
    finalize();

    return *this;
  }

  Index size(const Index outerIdx) const { return outer[outerIdx].size(); }

  Index nonZeros() const {
    Index nnz = 0;
    for (unsigned i = 0; i < outer.size(); ++i) nnz += outer[i].size();

    return nnz;
  }

  void changePtr(const InnerIterator& it, BlockPtr ptr) {
    const_cast<BlockPtr&>(it.asStdIterator()->second) = ptr;
  }
};

template <bool Compressed, typename Index_, typename BlockPtr_,
          template <typename> class ArrayType>
struct SparseBlockIndexTraits<
    SparseBlockIndex<Compressed, Index_, BlockPtr_, ArrayType> > {
  typedef Index_ Index;
  typedef BlockPtr_ BlockPtr;

  typedef SparseBlockIndex<Compressed, Index_, BlockPtr_, ArrayType>
      SparseBlockIndexType;

  //! Forward iterator
  struct InnerIterator {
    // Warning: This class does not implement the full RandomAccessIterator
    // concept ; only the operations that are required by std::lower_bound are
    // implemented
    typedef std::random_access_iterator_tag iterator_category;
    typedef Index value_type;
    typedef std::ptrdiff_t difference_type;
    typedef const Index* pointer;
    typedef const Index& reference;

    InnerIterator() {}

    InnerIterator(const SparseBlockIndexType& index, Index outer)
        : m_it(index.outer[outer].begin()), m_end(index.outer[outer].end()) {}

    operator bool() const { return m_it != m_end; }

    InnerIterator& operator++() {
      ++m_it;
      return *this;
    }
    InnerIterator& operator--() {
      --m_it;
      return *this;
    }

    InnerIterator& operator+=(const std::size_t n) {
      m_it += n;
      return *this;
    }

    difference_type operator-(const InnerIterator& other) const {
      return m_it - other.m_it;
    }

    Index operator*() const { return inner(); }

    InnerIterator end() const { return InnerIterator(*this).toEnd(); }

    InnerIterator& toEnd() {
      m_it = m_end;
      return *this;
    }

    bool after(Index outer) const { return inner() > outer; }

    Index inner() const { return m_it->first; }
    BlockPtr ptr() const { return m_it->second; }

    typename SparseBlockIndexType::Inner::const_iterator asStdIterator() const {
      return m_it;
    }

   private:
    typename SparseBlockIndexType::Inner::const_iterator m_it;
    typename SparseBlockIndexType::Inner::const_iterator m_end;
  };
};

}  // namespace bogus

#endif  // SPARSEBLOCKINDEX_HPP
