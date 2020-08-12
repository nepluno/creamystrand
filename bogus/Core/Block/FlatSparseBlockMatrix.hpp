/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_FLAT_SPARSEBLOCKMATRIX_HPP
#define BOGUS_FLAT_SPARSEBLOCKMATRIX_HPP

#include "SparseBlockMatrixBase.hpp"

namespace bogus {

namespace internal {

//! Storage for blocks shapes and coefficient data of FlatSparseBlockMatrix
template <typename BlockT, typename Index, typename BlockPtr>
struct FlatBlockStorage {
  typedef BlockTraits<BlockT> Traits;
  typedef typename Traits::Scalar Scalar;

  typedef typename BlockT::MapType BlockRef;
  typedef typename BlockT::ConstMapType ConstBlockRef;

  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Buffer;

  //! Shape and offset in data buffer for a given block
  struct Block {
    BlockPtr offset;
    Index rows;
    Index cols;

    explicit Block(BlockPtr o = -1, Index r = Traits::RowsAtCompileTime,
                   Index c = Traits::ColsAtCompileTime)
        : offset(o), rows(r), cols(c) {}

    void resize(Index r, Index c) {
      rows = r;
      cols = c;
    }

#ifdef BOGUS_WITH_BOOST_SERIALIZATION
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar& offset;
      ar& rows;
      ar& cols;
    }
#endif
  };

  FlatBlockStorage() : m_dataSize(0) {}

  // Buffer

  template <typename BlocksType>
  static BlockPtr totalNonZeros(const BlocksType& blocks) {
    BlockPtr totSize = 0;
    for (BlockPtr i = 0; i < blocks.size(); ++i) {
      totSize += blocks[i].rows() * blocks[i].cols();
    }
    return totSize;
  }

  void resetBuffer(std::size_t size) { m_buffer.resize(size); }

  void resizeBuffer(std::size_t size) {
    Buffer newBuf(size);
    newBuf.head(std::min(m_buffer.size(), newBuf.size())) = m_buffer;
    m_buffer.swap(newBuf);
  }

  void grow(Index n) {
    m_dataSize += n;
    if (m_dataSize >= m_buffer.size()) resizeBuffer((3 * m_dataSize) / 2);
  }

  void fitBufferToDataSize() { m_buffer.resize(m_dataSize); }

  // Block shapes

  void reserve(const BlockPtr nBlocks) { m_shapes.reserve(nBlocks); }

  BlockPtr append(const Index rows, const Index cols) {
    BlockPtr ptr = size();
    m_shapes.push_back(Block(m_dataSize, rows, cols));

    grow(rows * cols);
    return ptr;
  }

  void preallocShapes(const Index n) {
    clear();
    m_shapes.resize(n);
  }

  //! Sets the rows and cols of a preallocated shape, increase tot NNZ
  void makeShape(const BlockPtr ptr, const Index rows, const Index cols) {
    m_shapes[ptr] = Block(m_dataSize, rows, cols);
    m_dataSize += rows * cols;
  }

  //! Copy shapes from block list and reset buffer
  template <bool Transpose, typename BlocksType>
  void fit(const BlocksType& blocks) {
    clear();

    for (BlockPtr i = 0; i < blocks.size(); ++i) {
      if (Transpose)
        append(blocks[i].cols(), blocks[i].rows());
      else
        append(blocks[i].rows(), blocks[i].cols());
    }

    resetBuffer(totalNonZeros());
  }

  // Cleanup

  void clear() {
    m_shapes.clear();
    m_dataSize = 0;
  }

  void swap(FlatBlockStorage& o) {
    using std::swap;
    swap(m_dataSize, o.m_dataSize);
    swap(m_shapes, o.m_shapes);
    m_buffer.swap(o.m_buffer);
  }

  // Acces
  BlockRef operator[](BlockPtr ptr) {
    const Block& block = m_shapes[ptr];
    return BlockRef(m_buffer.data() + block.offset, block.rows, block.cols);
  }
  ConstBlockRef operator[](BlockPtr ptr) const {
    const Block& block = m_shapes[ptr];
    return ConstBlockRef(m_buffer.data() + block.offset, block.rows,
                         block.cols);
  }

  BlockPtr size() const { return m_shapes.size(); }

  typename Buffer::Index totalNonZeros() const { return m_dataSize; }

  const Buffer& buffer() const { return m_buffer; }

  typedef typename ResizableSequenceContainer<Block>::Type BlocksArrayType;
  const BlocksArrayType& shapes() const { return m_shapes; }

#ifdef BOGUS_WITH_BOOST_SERIALIZATION
  template <typename Archive>
  void serialize(Archive& ar, const unsigned int) {
    ar& m_buffer;
    ar& m_dataSize;
    ar& m_shapes;
  }
#endif

 protected:
  Buffer m_buffer;
  typename Buffer::Index m_dataSize;

  BlocksArrayType m_shapes;
};

}  // namespace internal

//! Specialization of BlockMatrixTraits for SparseBlockMatrix
template <typename BlockT, int Flags>
struct BlockMatrixTraits<FlatSparseBlockMatrix<BlockT, Flags> >
    : public BlockMatrixTraits<
          BlockObjectBase<FlatSparseBlockMatrix<BlockT, Flags> > > {
  typedef BlockMatrixTraits<
      BlockObjectBase<FlatSparseBlockMatrix<BlockT, Flags> > >
      BaseTraits;
  typedef typename BaseTraits::Index Index;
  typedef typename BlockTraits<BlockT>::Scalar Scalar;

  enum {
    is_symmetric = !!(Flags & flags::SYMMETRIC),
    RowsPerBlock = BlockTraits<BlockT>::RowsAtCompileTime,
    ColsPerBlock = BlockTraits<BlockT>::ColsAtCompileTime
  };

  typedef BOGUS_DEFAULT_BLOCK_PTR_TYPE BlockPtr;
  typedef BlockT BlockType;

  template <typename OtherBlockType>
  struct ResizableBlockContainer {
    typedef internal::FlatBlockStorage<OtherBlockType, Index, BlockPtr> Type;
  };

  typedef typename ResizableBlockContainer<BlockType>::Type BlocksArrayType;

  typedef typename BlocksArrayType::BlockRef BlockRef;
  typedef typename BlocksArrayType::ConstBlockRef ConstBlockRef;

  enum {
    is_compressed = !!(~Flags & flags::UNCOMPRESSED),
    is_col_major = !!(Flags & flags::COL_MAJOR),
    flags = Flags
  };

  typedef SparseBlockIndex<is_compressed, Index, BlockPtr> MajorIndexType;
};

//! Sparse Block Matrix
/*!
  \tparam BlockT the type of the blocks of the matrix. Can be scalar, Eigen
  dense of sparse matrices, or basically anything provided a few functions are
  specialized \tparam Flags a combination of the values defined in \ref
  bogus::flags
  */
template <typename BlockT, int Flags>
class FlatSparseBlockMatrix
    : public SparseBlockMatrixBase<FlatSparseBlockMatrix<BlockT, Flags> > {
 public:
  typedef SparseBlockMatrixBase<FlatSparseBlockMatrix> Base;
  typedef typename Base::Scalar Scalar;

  FlatSparseBlockMatrix() : Base() {}

  template <typename Index>
  FlatSparseBlockMatrix(Index rowsOfBlocks, Index colsOfBlocks) : Base() {
    BOGUS_STATIC_ASSERT(Base::has_fixed_size_blocks,
                        BLOCKS_MUST_HAVE_FIXED_DIMENSIONS);
    Base::setRows(rowsOfBlocks);
    Base::setCols(colsOfBlocks);
  }

  template <typename RhsT>
  FlatSparseBlockMatrix(const BlockObjectBase<RhsT>& rhs) : Base() {
    Base::operator=(rhs.derived());
  }

  template <typename RhsT>
  FlatSparseBlockMatrix& operator=(const BlockObjectBase<RhsT>& rhs) {
    return (Base::operator=(rhs.derived())).derived();
  }

  // Allocations and stuff

  using Base::m_blocks;
  using Base::m_transposeBlocks;
  typedef typename Base::Traits Traits;
  typedef typename Base::Index Index;
  typedef typename Base::MajorIndexType MajorIndexType;
  typedef typename Base::BlockPtr BlockPtr;

  template <bool EnforceThreadSafety>
  typename Base::BlockRef allocateBlock(typename Base::BlockPtr& ptr,
                                        Index outerSize, Index innerSize) {
#ifndef BOGUS_DONT_PARALLELIZE
    Lock::Guard<EnforceThreadSafety> guard(Base::m_lock);
#endif
    const Index cols = Traits::is_col_major ? outerSize : innerSize;
    const Index rows = Traits::is_col_major ? innerSize : outerSize;
    ptr = m_blocks.append(rows, cols);
    return m_blocks[ptr];
  }

  void reserve(std::size_t nBlocks, std::size_t totalNonZeros) {
    m_blocks.reserve(nBlocks);
    Base::majorIndex().reserve(nBlocks);
    m_blocks.resizeBuffer(totalNonZeros);
  }

  template <typename BlocksType>
  void resetFor(const BlocksType& blocks) {
    Base::clear();
    m_blocks.reserve(blocks.size());
    Base::majorIndex().reserve(blocks.size());
    m_blocks.resetBuffer(m_blocks.totalNonZeros(blocks));
  }

  template <bool Transpose, typename BlocksType>
  void copyBlockShapes(const BlocksType& blocks) {
    m_blocks.template fit<Transpose>(blocks);
  }

  template <bool Transpose, typename BlocksType>
  void copyTransposeBlockShapes(const BlocksType& blocks) {
    m_transposeBlocks.template fit<Transpose>(blocks);
  }

  template <typename IndexType>
  void createBlockShapes(const BlockPtr nBlocks, const IndexType& index,
                         typename Traits::BlocksArrayType& blocks) {
    blocks.preallocShapes(nBlocks);

    for (Index outer = 0; outer < index.outerSize(); ++outer) {
      for (typename IndexType::InnerIterator it(index, outer); it; ++it) {
        const Index inner = it.inner();
        const Index cols = Traits::is_col_major ? Base::blockCols(outer)
                                                : Base::blockCols(inner);
        const Index rows = Traits::is_col_major ? Base::blockRows(inner)
                                                : Base::blockRows(outer);

        blocks.makeShape(it.ptr(), rows, cols);
      }
    }
    blocks.resetBuffer(blocks.totalNonZeros());
  }
};

// Specialization for block matrix of MappedSparseBlockMatrix
template <typename BlockT, int Flags>
struct BlockTraits<FlatSparseBlockMatrix<BlockT, Flags> > {
  typedef FlatSparseBlockMatrix<BlockT, Flags> BlockType;

  typedef typename BlockType::Scalar Scalar;
  typedef FlatSparseBlockMatrix<typename BlockType::TransposeBlockType,
                                Flags ^ flags::COL_MAJOR>
      TransposeStorageType;

  enum {
    RowsAtCompileTime = internal::DYNAMIC,
    ColsAtCompileTime = internal::DYNAMIC,
    uses_plain_array_storage = 0,
    is_row_major = !BlockMatrixTraits<BlockType>::is_col_major,
    is_self_transpose = BlockMatrixTraits<BlockType>::is_symmetric
  };
};

}  // namespace bogus

#endif  // SPARSEBLOCKMATRIX_HH
