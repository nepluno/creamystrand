/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*! \file
        Necessary bindings to use Eigen matrices as block types, and
        \c operator* specialization for matrix/vector products
*/

#ifndef BLOCK_EIGENBINDINGS_HPP
#define BLOCK_EIGENBINDINGS_HPP

#include <Eigen/Core>

#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE
#include "SparseHeader.hpp"
#endif

#include "../Block/BlockMatrixBase.hpp"
#include "../Block/Expressions.hpp"

#ifndef BOGUS_BLOCK_WITHOUT_LINEAR_SOLVERS
#include "EigenLinearSolvers.hpp"
#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE
#include "EigenSparseLinearSolvers.hpp"
#endif
#endif

#include "../Utils/CppTools.hpp"

#define BOGUS_EIGEN_NEW_EXPRESSIONS EIGEN_VERSION_AT_LEAST(3, 2, 90)

namespace bogus {

// transpose_block, is_zero, resize, set_identity

template <typename EigenDerived>
inline bool is_zero(const Eigen::MatrixBase<EigenDerived>& block,
                    typename EigenDerived::Scalar precision) {
  return block.isZero(precision);
}

template <typename EigenDerived>
inline void set_zero(Eigen::MatrixBase<EigenDerived>& block) {
  block.derived().setZero();
}

template <typename EigenDerived>
inline void set_identity(Eigen::MatrixBase<EigenDerived>& block) {
  block.derived().setIdentity();
}

template <typename EigenDerived>
inline void resize(Eigen::MatrixBase<EigenDerived>& block, int rows, int cols) {
  block.derived().resize(rows, cols);
}

template <typename EigenDerived>
inline const typename EigenDerived::Scalar* data_pointer(
    const Eigen::MatrixBase<EigenDerived>& block) {
  return block.derived().data();
}

#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE

#if !BOGUS_EIGEN_NEW_EXPRESSIONS
template <typename BlockT>
struct BlockTransposeTraits<Eigen::SparseMatrixBase<BlockT> > {
  typedef const Eigen::Transpose<const BlockT> ReturnType;
};
template <typename _Scalar, int _Flags, typename _Index>
struct BlockTransposeTraits<Eigen::SparseMatrix<_Scalar, _Flags, _Index> >
    : public BlockTransposeTraits<Eigen::SparseMatrixBase<
          Eigen::SparseMatrix<_Scalar, _Flags, _Index> > > {};
template <typename _Scalar, int _Flags, typename _Index>
struct BlockTransposeTraits<Eigen::SparseVector<_Scalar, _Flags, _Index> >
    : public BlockTransposeTraits<Eigen::SparseMatrixBase<
          Eigen::SparseVector<_Scalar, _Flags, _Index> > > {};
template <typename _Scalar, int _Flags, typename _Index>
struct BlockTransposeTraits<Eigen::MappedSparseMatrix<_Scalar, _Flags, _Index> >
    : public BlockTransposeTraits<Eigen::SparseMatrixBase<
          Eigen::MappedSparseMatrix<_Scalar, _Flags, _Index> > > {};

template <typename EigenDerived>
inline const Eigen::Transpose<const EigenDerived> transpose_block(
    const Eigen::SparseMatrixBase<EigenDerived>& block) {
  return block.transpose();
}
#endif

template <typename EigenDerived>
inline bool is_zero(const Eigen::SparseMatrixBase<EigenDerived>& block,
                    typename EigenDerived::Scalar precision) {
  return block.isZero(precision);
}

template <typename Scalar, int Options, typename Index>
inline void set_identity(Eigen::SparseMatrix<Scalar, Options, Index>& block) {
  return block.setIdentity();
}

template <typename Scalar, int Options, typename Index>
inline void resize(Eigen::SparseMatrix<Scalar, Options, Index>& block,
                   Index rows, Index cols) {
  block.resize(rows, cols);
}

#endif

// Block traits for Eigen::Matrix

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows,
          int _MaxCols>
struct BlockTraits<
    Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> > {
  typedef _Scalar Scalar;
  enum {
    RowsAtCompileTime = _Rows,
    ColsAtCompileTime = _Cols,
    uses_plain_array_storage = 1,
    is_row_major = !!(_Options & Eigen::RowMajor),
    is_self_transpose = (_Rows == _Cols) && (_Rows == 1)
  };

  // Manipulates _Options so that row and column vectorz have the correct
  // RowMajor value Should we default to _Options ^ Eigen::RowMajor ?
  typedef Eigen::Matrix<_Scalar, _Cols, _Rows,
                        (_Options |
                         ((_Cols == 1 && _Rows != 1) ? Eigen::RowMajor : 0)) &
                            ~((_Rows == 1 && _Cols != 1) ? Eigen::RowMajor : 0),
                        _MaxCols, _MaxRows>
      TransposeStorageType;
};

// Block/block product return type

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows,
          int _MaxCols, typename _Scalar2, int _Rows2, int _Cols2,
          int _Options2, int _MaxRows2, int _MaxCols2, bool TransposeLhs,
          bool TransposeRhs>
struct BlockBlockProductTraits<
    Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>,
    Eigen::Matrix<_Scalar2, _Rows2, _Cols2, _Options2, _MaxRows2, _MaxCols2>,
    TransposeLhs, TransposeRhs> {
  typedef Eigen::Matrix<_Scalar, SwapIf<TransposeLhs, _Rows, _Cols>::First,
                        SwapIf<TransposeRhs, _Rows2, _Cols2>::Second, _Options,
                        SwapIf<TransposeLhs, _MaxRows, _MaxCols>::First,
                        SwapIf<TransposeRhs, _MaxRows2, _MaxCols2>::Second>
      ReturnType;
};

template <typename Derived>
inline typename Eigen::internal::plain_matrix_type<Derived>::type
get_mutable_vector(const Eigen::MatrixBase<Derived>&) {
  return typename Eigen::internal::plain_matrix_type<Derived>::type();
}

}  // namespace bogus

// Matrix vector product return types and operator*

namespace bogus {
namespace mv_impl {

//! Wrapper so our SparseBlockMatrix can be used inside Eigen expressions
template <typename Derived>
struct EigenBlockWrapper
    : public Eigen::EigenBase<EigenBlockWrapper<Derived> > {
  typedef EigenBlockWrapper Nested;
  typedef EigenBlockWrapper NestedExpression;
  typedef EigenBlockWrapper PlainObject;
  typedef Eigen::internal::traits<EigenBlockWrapper<Derived> > Traits;
  typedef typename Traits::Scalar Scalar;
  typedef typename Traits::Index Index;

  enum {
    Flags = Traits::Flags,
    IsVectorAtCompileTime = 0,
    MaxRowsAtCompileTime = Traits::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = Traits::MaxColsAtCompileTime
  };

  EigenBlockWrapper(const Derived& obj_, Scalar s = 1)
      : obj(obj_), scaling(s) {}

  Index rows() const { return obj.rows(); }
  Index cols() const { return obj.cols(); }

  // Mult by scalar

  inline EigenBlockWrapper operator*(const Scalar& scalar) const {
    return EigenBlockWrapper(obj, scaling * scalar);
  }

  inline friend EigenBlockWrapper operator*(const Scalar& scalar,
                                            const EigenBlockWrapper& matrix) {
    return EigenBlockWrapper(matrix.obj, matrix.scaling * scalar);
  }

  // Product with other Eigen expr
#if BOGUS_EIGEN_NEW_EXPRESSIONS
  template <typename EigenDerived>
  inline Eigen::Product<EigenBlockWrapper, EigenDerived> operator*(
      const EigenDerived& rhs) const {
    return Eigen::Product<EigenBlockWrapper, EigenDerived>(*this, rhs);
  }
  template <typename EigenDerived>
  friend inline Eigen::Product<EigenDerived, EigenBlockWrapper> operator*(
      const EigenDerived& lhs, const EigenBlockWrapper& matrix) {
    return Eigen::Product<EigenDerived, EigenBlockWrapper>(lhs, matrix);
  }
#endif

  const Derived& obj;
  const Scalar scaling;
};

//! SparseBlockMatrix / Dense producty expression

template <typename Lhs, typename Rhs>
struct block_product_impl;

#if BOGUS_EIGEN_NEW_EXPRESSIONS
template <typename Lhs, typename Rhs>
struct BlockEigenProduct : public Eigen::Product<Lhs, Rhs> {
  typedef Eigen::Product<Lhs, Rhs> Base;

  BlockEigenProduct(const Lhs& lhs, const Rhs& rhs) : Base(lhs, rhs) {}
};
#else

template <typename Lhs, typename Rhs>
struct BlockEigenProduct
    : public Eigen::ProductBase<BlockEigenProduct<Lhs, Rhs>, Lhs, Rhs> {
  typedef Eigen::ProductBase<BlockEigenProduct, Lhs, Rhs> Base;

  BlockEigenProduct(const Lhs& lhs, const Rhs& rhs) : Base(lhs, rhs) {}

  EIGEN_DENSE_PUBLIC_INTERFACE(BlockEigenProduct)
  using Base::m_lhs;
  using Base::m_rhs;
  typedef block_product_impl<typename Base::LhsNested, typename Base::RhsNested>
      product_impl;

  template <typename Dest>
  inline void evalTo(Dest& dst) const {
    product_impl::evalTo(dst, m_lhs, this->m_rhs);
  }
  template <typename Dest>
  inline void scaleAndAddTo(Dest& dst, Scalar alpha) const {
    product_impl::scaleAndAddTo(dst, m_lhs, m_rhs, alpha);
  }
};
#endif

template <typename Derived, typename EigenDerived>
BlockEigenProduct<EigenBlockWrapper<Derived>, EigenDerived> block_eigen_product(
    const Derived& matrix, const EigenDerived& vector,
    typename Derived::Scalar scaling = 1) {
  return BlockEigenProduct<EigenBlockWrapper<Derived>, EigenDerived>(
      EigenBlockWrapper<Derived>(matrix, scaling), vector);
}

template <typename Derived, typename EigenDerived>
BlockEigenProduct<EigenDerived, EigenBlockWrapper<Derived> >
eigen_block_product(const EigenDerived& vector, const Derived& matrix,
                    typename Derived::Scalar scaling = 1) {
  return BlockEigenProduct<EigenDerived, EigenBlockWrapper<Derived> >(
      vector, EigenBlockWrapper<Derived>(matrix, scaling));
}

template <typename Derived, typename EigenDerived>
struct block_product_impl<bogus::mv_impl::EigenBlockWrapper<Derived>,
                          EigenDerived> {
  typedef bogus::mv_impl::EigenBlockWrapper<Derived> Lhs;
  typedef EigenDerived Rhs;
  typedef typename Derived::Scalar Scalar;

  template <typename Dst>
  static void evalTo(Dst& dst, const Lhs& lhs, const Rhs& rhs) {
    lhs.obj.template multiply<false>(rhs.derived(), dst, lhs.scaling, 0);
  }
  template <typename Dst>
  static void scaleAndAddTo(Dst& dst, const Lhs& lhs, const Rhs& rhs,
                            const Scalar& alpha) {
    lhs.obj.template multiply<false>(rhs.derived(), dst, alpha * lhs.scaling,
                                     1);
  }
};

template <typename Derived, typename EigenDerived>
struct block_product_impl<EigenDerived,
                          bogus::mv_impl::EigenBlockWrapper<Derived> > {
  typedef EigenDerived Lhs;
  typedef bogus::mv_impl::EigenBlockWrapper<Derived> Rhs;
  typedef typename Derived::Scalar Scalar;

  template <typename Dst>
  static void evalTo(Dst& dst, const Lhs& lhs, const Rhs& rhs) {
    Eigen::Transpose<Dst> transposed(dst.transpose());
    rhs.obj.template multiply<true>(lhs.transpose(), transposed, rhs.scaling,
                                    0);
  }
  template <typename Dst>
  static void scaleAndAddTo(Dst& dst, const Lhs& lhs, const Rhs& rhs,
                            const Scalar& alpha) {
    Eigen::Transpose<Dst> transposed(dst.transpose());
    rhs.obj.template multiply<true>(lhs.transpose(), transposed,
                                    alpha * rhs.scaling, 1);
  }
};

}  // namespace mv_impl
}  // namespace bogus

#if BOGUS_EIGEN_NEW_EXPRESSIONS
namespace Eigen {
namespace internal {

template <typename Derived, typename EigenDerived, int ProductType>
struct generic_product_impl<bogus::mv_impl::EigenBlockWrapper<Derived>,
                            EigenDerived, SparseShape, DenseShape, ProductType>
    : public generic_product_impl_base<
          bogus::mv_impl::EigenBlockWrapper<Derived>, EigenDerived,
          generic_product_impl<bogus::mv_impl::EigenBlockWrapper<Derived>,
                               EigenDerived, SparseShape, DenseShape,
                               ProductType> > {
  typedef bogus::mv_impl::EigenBlockWrapper<Derived> Lhs;
  typedef EigenDerived Rhs;
  typedef typename Derived::Scalar Scalar;

  typedef typename nested_eval<Rhs, Dynamic>::type RhsNested;
  typedef bogus::mv_impl::block_product_impl<Lhs, RhsNested> product_impl;

  template <typename Dst>
  static void evalTo(Dst& dst, const Lhs& lhs, const Rhs& rhs) {
    RhsNested rhsNested(rhs);
    product_impl::evalTo(dst, lhs, rhsNested);
  }
  template <typename Dst>
  static void scaleAndAddTo(Dst& dst, const Lhs& lhs, const Rhs& rhs,
                            Scalar alpha) {
    RhsNested rhsNested(rhs);
    product_impl::scaleAndAddTo(dst, lhs, rhsNested, alpha);
  }
};

template <typename Derived, typename EigenDerived, int ProductType>
struct generic_product_impl<EigenDerived,
                            bogus::mv_impl::EigenBlockWrapper<Derived>,
                            DenseShape, SparseShape, ProductType>
    : public generic_product_impl_base<
          EigenDerived, bogus::mv_impl::EigenBlockWrapper<Derived>,
          generic_product_impl<EigenDerived,
                               bogus::mv_impl::EigenBlockWrapper<Derived>,
                               DenseShape, SparseShape, ProductType> > {
  typedef EigenDerived Lhs;
  typedef bogus::mv_impl::EigenBlockWrapper<Derived> Rhs;
  typedef typename Derived::Scalar Scalar;

  typedef typename nested_eval<Lhs, Dynamic>::type LhsNested;
  typedef bogus::mv_impl::block_product_impl<LhsNested, Rhs> product_impl;

  template <typename Dst>
  static void evalTo(Dst& dst, const Lhs& lhs, const Rhs& rhs) {
    LhsNested lhsNested(lhs);
    product_impl::evalTo(dst, lhsNested, rhs);
  }
  template <typename Dst>
  static void scaleAndAddTo(Dst& dst, const Lhs& lhs, const Rhs& rhs,
                            const Scalar& alpha) {
    LhsNested lhsNested(lhs);
    product_impl::scaleAndAddTo(dst, lhsNested, rhs, alpha);
  }
};

// s * (A * V) -> ( (s*A) * V )
// (A already includes a scaling parameter)
// TODO adapt to other orderings
template <typename Derived, typename Rhs, typename Scalar1, typename Scalar2,
          typename Plain1>
struct evaluator<CwiseBinaryOp<
    internal::scalar_product_op<Scalar1, Scalar2>,
    const CwiseNullaryOp<internal::scalar_constant_op<Scalar1>, Plain1>,
    const Product<bogus::mv_impl::EigenBlockWrapper<Derived>, Rhs,
                  DefaultProduct> > >
    : public evaluator<Product<bogus::mv_impl::EigenBlockWrapper<Derived>, Rhs,
                               DefaultProduct> > {
  typedef CwiseBinaryOp<
      internal::scalar_product_op<Scalar1, Scalar2>,
      const CwiseNullaryOp<internal::scalar_constant_op<Scalar1>, Plain1>,
      const Product<bogus::mv_impl::EigenBlockWrapper<Derived>, Rhs,
                    DefaultProduct> >
      XprType;
  typedef evaluator<
      Product<bogus::mv_impl::EigenBlockWrapper<Derived>, Rhs, DefaultProduct> >
      Base;

  explicit evaluator(const XprType& xpr)
      : Base(bogus::mv_impl::EigenBlockWrapper<Derived>(
                 xpr.rhs().lhs().obj,
                 xpr.lhs().functor().m_other * xpr.rhs().lhs().scaling) *
             xpr.rhs().rhs()) {}
};

}  // namespace internal
}  // namespace Eigen

#endif

// Eigen traits for our new structs
namespace Eigen {
namespace internal {

template <typename Lhs, typename Rhs>
struct traits<bogus::mv_impl::BlockEigenProduct<Lhs, Rhs> >
#if BOGUS_EIGEN_NEW_EXPRESSIONS
    : public traits<typename bogus::mv_impl::BlockEigenProduct<Lhs, Rhs>::Base>
#else
    : public traits<
          ProductBase<bogus::mv_impl::BlockEigenProduct<Lhs, Rhs>, Lhs, Rhs> >
#endif
{
  typedef Dense StorageKind;
};

template <typename Derived>
struct traits<bogus::mv_impl::EigenBlockWrapper<Derived> > {
  typedef typename Derived::Scalar Scalar;
  typedef typename Derived::Index Index;
  typedef typename Derived::Index StorageIndex;
#if BOGUS_EIGEN_NEW_EXPRESSIONS
  typedef Sparse StorageKind;
#else
  typedef Dense StorageKind;
#endif
  typedef MatrixXpr XprKind;
  enum {
    RowsAtCompileTime = Dynamic,
    ColsAtCompileTime = Dynamic,
    MaxRowsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic,
    Flags = 0
  };
};
}  // namespace internal

}  // namespace Eigen

// Matrix-Vector operators

namespace bogus {

template <typename Derived, typename EigenDerived>
mv_impl::BlockEigenProduct<mv_impl::EigenBlockWrapper<Derived>, EigenDerived>
operator*(const BlockObjectBase<Derived>& lhs,
          const Eigen::MatrixBase<EigenDerived>& rhs) {
  assert(rhs.rows() == lhs.cols());
  return mv_impl::block_eigen_product(lhs.derived(), rhs.derived());
}

template <typename Derived, typename EigenDerived>
mv_impl::BlockEigenProduct<mv_impl::EigenBlockWrapper<Derived>, EigenDerived>
operator*(const Scaling<Derived>& lhs,
          const Eigen::MatrixBase<EigenDerived>& rhs) {
  assert(rhs.rows() == lhs.cols());
  return mv_impl::block_eigen_product(lhs.operand.object, rhs.derived(),
                                      lhs.operand.scaling);
}

template <typename Derived, typename EigenDerived>
mv_impl::BlockEigenProduct<EigenDerived, mv_impl::EigenBlockWrapper<Derived> >
operator*(const Eigen::MatrixBase<EigenDerived>& lhs,
          const BlockObjectBase<Derived>& rhs) {
  assert(lhs.cols() == rhs.rows());
  return mv_impl::eigen_block_product(lhs.derived(), rhs.derived());
}

template <typename Derived, typename EigenDerived>
mv_impl::BlockEigenProduct<EigenDerived, mv_impl::EigenBlockWrapper<Derived> >
operator*(const Eigen::MatrixBase<EigenDerived>& lhs,
          const Scaling<Derived>& rhs) {
  assert(lhs.cols() == rhs.rows());
  return mv_impl::eigen_block_product(lhs.derived(), rhs.operand.object,
                                      rhs.operand.scaling);
}

}  // namespace bogus

#endif  // EIGENBINDINGS_HPP
