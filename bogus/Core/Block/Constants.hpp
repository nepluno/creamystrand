/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_CONSTANTS_HPP
#define BOGUS_BLOCK_CONSTANTS_HPP

#include <cstddef>

#if !(defined(_OPENMP) || defined(BOGUS_DONT_PARALLELIZE))
#define BOGUS_DONT_PARALLELIZE
#endif

#ifndef BOGUS_DEFAULT_INDEX_TYPE
#define BOGUS_DEFAULT_INDEX_TYPE int
#endif

#ifndef BOGUS_DEFAULT_BLOCK_PTR_TYPE
#define BOGUS_DEFAULT_BLOCK_PTR_TYPE std::size_t
#endif

#ifndef BOGUS_DEFAULT_DENSE_INDEX_TYPE
#ifndef EIGEN_DEFAULT_DENSE_INDEX_TYPE
#define BOGUS_DEFAULT_DENSE_INDEX_TYPE int
#else
#define BOGUS_DEFAULT_DENSE_INDEX_TYPE EIGEN_DEFAULT_DENSE_INDEX_TYPE
#endif
#endif

namespace bogus {

//! Flags for compile-time tuning of the behavior of objects such as
//! SparseBlockMatrix
/*! Any combination if those is possible, using the 'binary or' ('|') operator.
 */
namespace flags {
enum {
  //! Default value: the matrix will use an uncompressed index, will be
  //! row-major, and not symmetric
  NONE = 0,
  //! Use an uncompressed index
  /*! This removes some restrictions on the order in which elements can be
     inserted, but can be less efficient and will dsallow interoperability with
     other formats such as MKL's BSR. <b>If the matrix can be created in a
     compressed way, that is with all its elements inserted in order, you
     probably should not set the UNCOMPRESSED flag. </b> \sa
     SparseBlockMatrixBase::insert() and SparseBlockMatrixBase::insertBack()
  */
  UNCOMPRESSED = 0x1,
  //! Store and index blocks in a column major way
  COL_MAJOR = 0x2,
  ROW_MAJOR = 0,  //!< Alias for convenience
                  //! Store only half the matrix, or rather the triangular part
                  //! for which \c inner \c <= \c outer,
  /*! \c outer being the row and \c inner the column for row-major matrices.
          Linear algebra operations such as matrix-vector and matrix-matrix
     multiplication will work just like if the matrix was fully populated, but
     at a lower memory access cost
  */
  SYMMETRIC = 0x4
};
}
// Reduce verbosity of public API
using namespace flags;

namespace internal {
//! Indicates a value that is not known at compile time
enum { DYNAMIC = -1 };
}  // namespace internal

typedef BOGUS_DEFAULT_DENSE_INDEX_TYPE DenseIndexType;

}  // namespace bogus

#endif
