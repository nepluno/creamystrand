/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2015 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_DYNAMIC_EXPRESSIONS_HPP
#define BOGUS_BLOCK_DYNAMIC_EXPRESSIONS_HPP

#include <list>

#include "DynamicExpressions.hpp"
#include "Expressions.hpp"

namespace bogus {

//! Sum of n similar expressions
template <typename Expression>
struct NarySum : public BlockObjectBase<NarySum<Expression> > {
  typedef BlockObjectBase<NarySum<Expression> > Base;
  typedef BlockMatrixTraits<NarySum<Expression> > Traits;

  typedef typename Traits::Index Index;
  typedef typename Base::Scalar Scalar;

  // Using list as we dont assume existence of c++11 vector::emplace_back
  typedef std::list<Scaling<Expression> > Sum;
  Sum members;

  NarySum(const Index rows, const Index cols) : m_rows(rows), m_cols(cols) {}
  NarySum(const Expression &expr) : m_rows(expr.rows()), m_cols(expr.cols()) {
    (*this) += expr;
  }
  NarySum(const Scaling<Expression> &expr)
      : m_rows(expr.rows()), m_cols(expr.cols()) {
    (*this) += expr;
  }

  NarySum &operator+=(const Scaling<Expression> &expr) {
    assert(expr.rows() == m_rows && expr.cols() == m_cols);
    members.push_back(expr);
    return *this;
  }

  NarySum &operator-=(const Scaling<Expression> &expr) {
    return (*this) +=
           Scaling<Expression>(expr.operand.object, -expr.operand.scaling);
  }

  NarySum &operator+=(const Expression &expr) {
    return (*this) += Scaling<Expression>(expr, 1);
  }

  NarySum &operator-=(const Expression &expr) {
    return (*this) += Scaling<Expression>(expr, -1);
  }

  NarySum &operator+=(const NarySum<Expression> &other) {
    assert(other.rows() == m_rows && other.cols() == m_cols);
    members.insert(members.begin(), other.members.begin(), other.members.end());
    return *this;
  }

  NarySum &operator-=(const NarySum<Expression> &other) {
    assert(other.rows() == m_rows && other.cols() == m_cols);
    for (typename Sum::const_iterator it = other.members.begin();
         it != other.members.end(); ++it) {
      (*this) -= *it;
    }
    return *this;
  }

  typename Base::ConstTransposeReturnType transpose() const {
    typename Base::ConstTransposeReturnType transposed_sum(m_cols, m_rows);
    for (typename Sum::const_iterator it = members.begin(); it != members.end();
         ++it) {
      transposed_sum += it->transpose();
    }

    return transposed_sum;
  }

  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT &rhs, ResT &res, Scalar alpha = 1,
                Scalar beta = 0) const;

  bool empty() const { return members.begin() == members.end(); }

  Index rows() const { return m_rows; }
  Index cols() const { return m_cols; }

  Index rowsOfBlocks() const {
    assert(!empty());
    return members.front().colsOfBlocks();
  }
  Index colsOfBlocks() const {
    assert(!empty());
    return members.front().rowsOfBlocks();
  }
  Index blockRows(Index row) const {
    assert(!empty());
    return members.front().blockCols(row);
  }
  Index blockCols(Index col) const {
    assert(!empty());
    return members.front().blockRows(col);
  }
  const Index *rowOffsets() const {
    assert(!empty());
    return members.front().colOffsets();
  }
  const Index *colOffsets() const {
    assert(!empty());
    return members.front().rowOffsets();
  }

 private:
  // Used for consistency checks only
  const Index m_rows;
  const Index m_cols;
};

template <typename Expression>
struct BlockMatrixTraits<NarySum<Expression> >
    : public BlockMatrixTraits<Addition<Expression, Expression> > {
  enum { is_temporary = 0 };

  typedef typename Expression::PlainObjectType PlainObjectType;
  typedef NarySum<typename Expression::TransposeObjectType>
      ConstTransposeReturnType;
  typedef ConstTransposeReturnType TransposeObjectType;
};

}  // namespace bogus

#endif
