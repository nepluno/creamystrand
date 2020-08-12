/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EIGENSERIALIZATION_HH_
#define EIGENSERIALIZATION_HH_

#include <boost/serialization/split_free.hpp>
#include <boost/serialization/split_member.hpp>

#include "../Core/Definitions.hh"

namespace boost {
namespace serialization {

template <class Archive, typename ScalarT, int Rows, int Cols>
void save(Archive& ar, const Eigen::Matrix<ScalarT, Rows, Cols>& mat,
          const unsigned int version) {
  long rows = mat.rows();
  long cols = mat.cols();

  ar << rows << cols;
  for (long i = 0; i < mat.rows(); ++i)
    for (long j = 0; j < mat.cols(); ++j) ar << mat(i, j);
}

template <class Archive, typename ScalarT, int Rows, int Cols>
void load(Archive& ar, Eigen::Matrix<ScalarT, Rows, Cols>& mat,
          const unsigned int version) {
  long rows, cols;
  ar >> rows >> cols;
  mat.resize(rows, cols);
  for (long i = 0; i < mat.rows(); ++i)
    for (long j = 0; j < mat.cols(); ++j) ar >> mat(i, j);
}

template <class Archive, typename ScalarT, int Rows, int Cols>
inline void serialize(Archive& ar, Eigen::Matrix<ScalarT, Rows, Cols>& mat,
                      const unsigned int version) {
  split_free(ar, mat, version);
}

template <class Archive, typename ScalarT>
void save(Archive& ar,
          const Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>& mat,
          const unsigned int version) {
  long rows = mat.rows();
  long cols = mat.cols();

  ar << rows << cols;
  for (long i = 0; i < mat.rows(); ++i)
    for (long j = 0; j < mat.cols(); ++j) ar << mat(i, j);
}

template <class Archive, typename ScalarT>
void load(Archive& ar,
          Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>& mat,
          const unsigned int version) {
  long rows, cols;
  ar >> rows >> cols;
  mat.resize(rows, cols);
  for (long i = 0; i < mat.rows(); ++i)
    for (long j = 0; j < mat.cols(); ++j) ar >> mat(i, j);
}

template <class Archive, typename ScalarT>
inline void serialize(
    Archive& ar, Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>& mat,
    const unsigned int version) {
  split_free(ar, mat, version);
}

template <class Archive, typename ScalarT>
void save(Archive& ar, const Eigen::Matrix<ScalarT, Eigen::Dynamic, 1>& vec,
          const unsigned int version) {
  const long size = vec.size();
  ar << size;
  for (long j = 0; j < size; ++j) ar << vec[j];
}

template <class Archive, typename ScalarT>
void load(Archive& ar, Eigen::Matrix<ScalarT, Eigen::Dynamic, 1>& vec,
          unsigned int version) {
  long size;
  ar >> size;
  vec.resize(size);
  for (long j = 0; j < size; ++j) ar >> vec[j];
}

template <class Archive, typename ScalarT>
inline void serialize(Archive& ar,
                      Eigen::Matrix<ScalarT, Eigen::Dynamic, 1>& vec,
                      const unsigned int version) {
  split_free(ar, vec, version);
}

template <class Archive, typename ScalarT, int m>
void save(Archive& ar, const Eigen::Matrix<ScalarT, m, 1>& vec,
          const unsigned int version) {
  for (long j = 0; j < m; ++j) ar << vec[j];
}

template <class Archive, typename ScalarT, int m>
void load(Archive& ar, Eigen::Matrix<ScalarT, m, 1>& vec,
          unsigned int version) {
  for (long j = 0; j < m; ++j) ar >> vec[j];
}

template <class Archive, typename ScalarT, int m>
inline void serialize(Archive& ar, Eigen::Matrix<ScalarT, m, 1>& vec,
                      const unsigned int version) {
  split_free(ar, vec, version);
}

template <class Archive, typename ScalarT, int _Options, typename _Index>
void save(Archive& ar,
          const Eigen::SparseMatrix<ScalarT, _Options, _Index>& mat,
          const unsigned int version) {
  long rows = mat.rows();
  long cols = mat.cols();
  long nnz = mat.nonZeros();

  ar << rows << cols << nnz;

  long inner;

  for (_Index outer = 0; outer < mat.outerSize(); ++outer) {
    nnz = mat.innerVector(outer).nonZeros();
    ar << nnz;

    typename Eigen::SparseMatrix<ScalarT, _Options, _Index>::InnerIterator it(
        mat, outer);
    for (; it; ++it) {
      inner = it.index();
      ar << inner << it.value();
    }
  }
}

template <class Archive, typename ScalarT, int _Options, typename _Index>
void load(Archive& ar, Eigen::SparseMatrix<ScalarT, _Options, _Index>& mat,
          const unsigned int version) {
  long rows, cols, nnz;
  ar >> rows >> cols >> nnz;
  mat.resize(rows, cols);
  mat.reserve(nnz);

  long inner, outer;
  ScalarT value;

  for (_Index outer = 0; outer < mat.outerSize(); ++outer) {
    ar >> nnz;
    for (long n = 0; n < nnz; ++n) {
      ar >> inner >> value;
      mat.insertByOuterInner(outer, inner) = value;
    }
  }
}

// template<class Archive, typename ScalarT, int _Options, typename _Index>
// inline void serialize( Archive & ar, Eigen::SparseMatrix<ScalarT, _Options,
// _Index> & mat
//         , const unsigned int version )
// {
//     split_free( ar, mat, version );
// }

}  // namespace serialization
}  // namespace boost

#endif /* EIGENSERIALIZATION_HH_ */
