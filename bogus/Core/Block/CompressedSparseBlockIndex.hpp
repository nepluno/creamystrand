/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef COMPRESSED_SPARSE_BLOCK_INDEX_HPP
#define COMPRESSED_SPARSE_BLOCK_INDEX_HPP

#include "../Utils/CppTools.hpp"
#include "SparseBlockIndex.hpp"

namespace bogus {

//! Compressed index, compatible with usual BSR/BSC formats
template <typename Index_, typename BlockPtr_,
          template <typename> class ArrayType>
struct SparseBlockIndex<true, Index_, BlockPtr_, ArrayType>
    : public SparseBlockIndexBase<
          SparseBlockIndex<true, Index_, BlockPtr_, ArrayType> > {
  typedef Index_ Index;
  typedef BlockPtr_ BlockPtr;

  typedef SparseBlockIndexBase<
      SparseBlockIndex<true, Index_, BlockPtr_, ArrayType> >
      Base;
  typedef typename Base::InnerOffsetsType InnerOffsetsType;
  typedef typename Base::InnerIterator InnerIterator;
  using Base::valid;

  typedef typename ArrayType<Index>::Type Inner;
  typedef typename ArrayType<Index>::Type Outer;

  InnerOffsetsType innerOffsets;

  //! Vector of inner indices
  Inner inner;
  //! Vector encoding the start and end of inner vectors
  Outer outer;

  SparseBlockIndex() : Base() {}

  void resizeOuter(Index size) { outer.assign(size + 1, 0); }
  void reserve(Index nnz) { inner.reserve(nnz); }

  Index outerSize() const { return outer.size() - 1; }
  Index nonZeros() const { return inner.size(); }

  const InnerOffsetsType& innerOffsetsArray() const { return innerOffsets; }

  //! \warning Only works for ordered insertion, and a call to \ref finalize()
  //! is always required once insertion is finished
  template <bool Ordered>
  void insert(Index outIdx, Index inIdx, BlockPtr ptr) {
    BOGUS_STATIC_ASSERT(Ordered, UNORDERED_INSERTION_WITH_COMPRESSED_INDEX);

    valid &= (ptr == (BlockPtr)(inner.size())) &&
             (0 == outer[outIdx + 1] || inIdx > inner.back());
    ++outer[outIdx + 1];
    inner.push_back(inIdx);
  }
  void insertBack(Index outIdx, Index inIdx, BlockPtr ptr) {
    insert<true>(outIdx, inIdx, ptr);
  }

  //! Finalizes the outer indices vector
  /*! Before calling this function, \c outer[i] contains the number of blocks
           in the \c i th inner vector

           After calling this functions, \c outer[i] contains the index of the
     start of the \c i th inner vector in \ref inner, and \c outer[i+1] its end
  */
  void finalize() {
    for (unsigned i = 1; i < outer.size(); ++i) {
      outer[i] += outer[i - 1];
    }
  }

  void clear() {
    outer.assign(outer.size(), 0);
    inner.clear();

    valid = true;
  }

  SparseBlockIndex& operator=(const SparseBlockIndex& compressed) {
    if (&compressed != this) {
      outer = compressed.outer;
      inner = compressed.inner;
      if (!compressed.innerOffsets.empty())
        innerOffsets = compressed.innerOffsets;
      valid = compressed.valid;
    }
    return *this;
  }

  template <typename SourceDerived>
  SparseBlockIndex& operator=(
      const SparseBlockIndexBase<SourceDerived>& source) {
    resizeOuter(source.outerSize());
    reserve(source.nonZeros());

    inner.clear();
    if (source.hasInnerOffsets()) {
      innerOffsets.resize(source.innerOffsetsArray().size());
      std::copy(source.innerOffsetsArray().begin(),
                source.innerOffsetsArray().end(), innerOffsets.begin());
    }
    valid = source.valid;

    for (typename SourceDerived::Index i = 0; i < source.outerSize(); ++i) {
      for (typename SourceDerived::InnerIterator it(source.derived(), i); it;
           ++it) {
        insertBack(i, it.inner(), it.ptr());
      }
    }

    finalize();

    return *this;
  }

  SparseBlockIndex& move(SparseBlockIndex& compressed) {
    if (&compressed != this) {
      // Want to have fun with gcc 4.6.3 ? Just swap the following statements !
      // Note: Swapping works perfectly with (at least) gcc 4.6.3 in
      // non-optimized mode, gcc 4.8, gcc 4.2, clang 4.2
      inner.swap(compressed.inner);
      outer.swap(compressed.outer);

      if (!compressed.innerOffsets.empty())
        innerOffsets.swap(compressed.innerOffsets);
      valid = compressed.valid;
      compressed.valid = false;
    }
    return *this;
  }

  template <typename SourceDerived>
  SparseBlockIndex& move(const SparseBlockIndexBase<SourceDerived>& source) {
    return (*this = source);
  }

  Index size(const Index outerIdx) const {
    return outer[outerIdx + 1] - outer[outerIdx];
  }

  // MKL BSR
  const Index* rowIndex() const { return data_pointer(outer); }
  const Index* columns() const { return data_pointer(inner); }

  // Same with nicer names
  const Index* outerIndexPtr() const { return data_pointer(outer); }
  const Index* innerIndexPtr() const { return data_pointer(inner); }
};

template <typename Index_, typename BlockPtr_,
          template <typename> class ArrayType>
struct SparseBlockIndexTraits<
    SparseBlockIndex<true, Index_, BlockPtr_, ArrayType> > {
  typedef Index_ Index;
  typedef BlockPtr_ BlockPtr;

  typedef SparseBlockIndex<true, Index_, BlockPtr_, ArrayType>
      SparseBlockIndexType;

  //! Forward iterator
  struct InnerIterator {
    // Warning: This class does not implement the full RandomAccessIterator
    // concept ; only the operations that are required by std::lower_bound are
    // implemented
    typedef std::random_access_iterator_tag iterator_category;
    typedef Index value_type;
    typedef ptrdiff_t difference_type;
    typedef const Index* pointer;
    typedef const Index& reference;

    InnerIterator() : m_inner(BOGUS_NULL_PTR(const Index)) {}

    InnerIterator(const SparseBlockIndexType& index, Index outer)
        : m_it(index.outer[outer]),
          m_end(index.outer[outer + 1]),
          m_inner(data_pointer(index.inner)) {}

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
      m_it = std::min(m_it + (Index)n, m_end);
      return *this;
    }

    difference_type operator-(const InnerIterator& other) const {
      return ((difference_type)m_it) - (difference_type)other.m_it;
    }

    Index operator*() const { return inner(); }

    InnerIterator end() const { return InnerIterator(*this).toEnd(); }

    InnerIterator& toEnd() {
      m_it = m_end;
      return *this;
    }

    bool after(Index outer) const { return inner() > outer; }

    Index inner() const { return m_inner[m_it]; }
    BlockPtr ptr() const { return m_it; }

    BlockPtr rawIndex() const { return m_it; }

   private:
    Index m_it;
    Index m_end;
    const Index* m_inner;
  };
};

}  // namespace bogus

#endif
