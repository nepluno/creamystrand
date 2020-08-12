/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef S_DEFINITIONS_HH
#define S_DEFINITIONS_HH

#ifdef WIN32
#define NOMINMAX
#endif

#include <stdint.h>

#include <iostream>
#include <memory>

#ifndef M_PI
#define M_PI 3.141592653589793238462650288
#endif
#ifndef M_PI_2
#define M_PI_2 (M_PI / 2.)
#endif
#ifndef M_PI_4
#define M_PI_4 (M_PI / 4.)
#endif

#ifndef STRANDSIM_INCLUDE_VANILLA_EIGEN
#define EIGEN_VECTOR_IO_FORMAT \
  Eigen::IOFormat(8, Eigen::DontAlignCols, ", ", ", ", "", "", "{ ", " }")
#define EIGEN_MATRIX_IO_FORMAT \
  Eigen::IOFormat(8, 0, ", ", "\n", "{ ", " }", "{ ", " }")
#undef EIGEN_DEFAULT_IO_FORMAT  // < To silence some warnings about redefining
#define EIGEN_DEFAULT_IO_FORMAT EIGEN_VECTOR_IO_FORMAT

#undef EIGEN_INITIALIZE_MATRICES_BY_ZERO  // < To silence some warnings about
                                          // redefining
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#endif          // STRANDSIM_INCLUDE_VANILLA_EIGEN
#undef Success  // Conflicts with Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/StdVector>

// template<typename Derived> std::ostream&
// operator<<( std::ostream& os, const Eigen::MatrixBase<Derived>& v )
//{
//    if ( v.cols() == 1 )
//        os << v.format( EIGEN_VECTOR_IO_FORMAT );
//    else
//        os << v.format( EIGEN_MATRIX_IO_FORMAT );
//
//    return os;
//}

namespace Eigen {
template <typename _Scalar, int _Options, typename _Index>
class SparseMatrix;
}

namespace strandsim {

typedef double Scalar;  ///< the scalar type
// typedef float Scalar; ///< the scalar type

typedef int IndexType;
typedef unsigned long long uint64;
typedef long long int64;

typedef Eigen::Matrix<Scalar, 2, 1> Vec2x;  ///< 2d scalar vector
typedef Eigen::Matrix<Scalar, 3, 1> Vec3x;  ///< 3d scalar vector
typedef Eigen::Matrix<Scalar, 4, 1> Vec4x;  ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 5, 1> Vec5x;  ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 6, 1> Vec6x;  ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 7, 1> Vec7x;  ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 8, 1> Vec8x;  ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 9, 1> Vec9x;  ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 11, 1>
    Vec11x;  ///< 11d scalar vector (stencil for local forces)
typedef Eigen::Matrix<Scalar, 27, 2>
    Mat27x2x;  ///< 27x2d scalar vector (stencil for local forces)
typedef Eigen::Matrix<Scalar, 27, 3>
    Mat27x3x;  ///< 27x3d scalar vector (stencil for local forces)
typedef Eigen::Matrix<Scalar, 27, 4>
    Mat27x4x;  ///< 27x4d scalar vector (stencil for local forces)
typedef Eigen::Matrix<Scalar, 27, 5>
    Mat27x5x;  ///< 27x5d scalar vector (stencil for local forces)
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
    VecXx;  ///< arbitrary dimension scalar vector
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, 1>
    VecXu;  ///< arbitrary dimension unsigned vector
typedef Eigen::Matrix<int, Eigen::Dynamic, 1>
    VecXi;  ///< arbitrary dimension unsigned vector
typedef Eigen::Matrix<unsigned char, Eigen::Dynamic, 1>
    VecXuc;  ///< arbitrary dimension unsigned char vector

typedef Eigen::Matrix<Scalar, 1, 2> Vec2xT;  ///< 2d scalar vector
typedef Eigen::Matrix<Scalar, 1, 3> Vec3xT;  ///< 3d scalar vector
typedef Eigen::Matrix<Scalar, 1, 4> Vec4xT;  ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 1, 5> Vec5xT;  ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 1, 6> Vec6xT;  ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 1, 7> Vec7xT;  ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 1, 8> Vec8xT;  ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 1, 9> Vec9xT;  ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 1, 11>
    Vec11xT;  ///< 11d scalar vector (stencil for local forces)

typedef Eigen::Matrix<Scalar, 1, Eigen::Dynamic>
    VecXxT;  ///< arbitrary dimension scalar vector
typedef Eigen::Matrix<unsigned, 1, Eigen::Dynamic>
    VecXuT;  ///< arbitrary dimension unsigned vector
typedef Eigen::Matrix<int, 1, Eigen::Dynamic>
    VecXiT;  ///< arbitrary dimension unsigned vector
typedef Eigen::Matrix<unsigned char, 1, Eigen::Dynamic>
    VecXucT;  ///< arbitrary dimension unsigned char vector

typedef std::vector<Vec2x, Eigen::aligned_allocator<Vec2x> >
    Vec2xArray;  ///< an array of 2d scalar vectors
// typedef std::vector<Vec3x, Eigen::aligned_allocator<Vec3x> > Vec3xArray; ///<
// an array of 3d scalar vectors
typedef std::vector<Vec3x, Eigen::aligned_allocator<Vec3x> > Vec3xArray;
typedef std::vector<Vec4x, Eigen::aligned_allocator<Vec4x> > Vec4xArray;
typedef std::vector<Vec11x, Eigen::aligned_allocator<Vec11x> >
    Vec11xArray;  ///< an array of 11d scalar vectors

typedef Eigen::Matrix<int, 2, 1> Vec2i;
typedef Eigen::Matrix<int, 3, 1> Vec3i;
typedef Eigen::Matrix<int, 4, 1> Vec4i;
typedef Eigen::Matrix<int, 5, 1> Vec5i;  ///< 4d scalar vector
typedef Eigen::Matrix<int, 6, 1> Vec6i;  ///< 4d scalar vector
typedef Eigen::Matrix<int, 7, 1> Vec7i;  ///< 4d scalar vector
typedef Eigen::Matrix<int, 8, 1> Vec8i;  ///< 4d scalar vector
typedef Eigen::Matrix<int, 9, 1> Vec9i;  ///< 4d scalar vector
typedef Eigen::Matrix<int, 27, 2>
    Mat27x2i;  ///< 27x2d int vector (stencil for local forces)
typedef Eigen::Matrix<int, 27, 3>
    Mat27x3i;  ///< 27x3d int vector (stencil for local forces)
typedef Eigen::Matrix<int, 27, 4>
    Mat27x4i;  ///< 27x4d int vector (stencil for local forces)
typedef Eigen::Matrix<int, 27, 5>
    Mat27x5i;  ///< 27x5d int vector (stencil for local forces)

typedef std::vector<Vec3i, Eigen::aligned_allocator<Vec3i> >
    Vec3iArray;  ///< an array of 3d float vectors
typedef std::vector<Vec4i, Eigen::aligned_allocator<Vec4i> >
    Vec4iArray;  ///< an array of 4d float vectors
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
    MatXi;  ///< arbitrary dimension scalar matrix

typedef Eigen::Matrix<short, 2, 1> Vec2s;
typedef Eigen::Matrix<short, 3, 1> Vec3s;
typedef Eigen::Matrix<short, 4, 1> Vec4s;
typedef Eigen::Matrix<short, 5, 1> Vec5s;  ///< 4d scalar vector
typedef Eigen::Matrix<short, 6, 1> Vec6s;  ///< 4d scalar vector
typedef Eigen::Matrix<short, 7, 1> Vec7s;  ///< 4d scalar vector
typedef Eigen::Matrix<short, 8, 1> Vec8s;  ///< 4d scalar vector
typedef Eigen::Matrix<short, 9, 1> Vec9s;  ///< 4d scalar vector
typedef Eigen::Matrix<short, 27, 2>
    Mat27x2s;  ///< 27x2d int vector (stencil for local forces)
typedef Eigen::Matrix<short, 27, 3>
    Mat27x3s;  ///< 27x3d int vector (stencil for local forces)
typedef Eigen::Matrix<short, 27, 4>
    Mat27x4s;  ///< 27x4d int vector (stencil for local forces)
typedef Eigen::Matrix<short, 27, 5>
    Mat27x5s;  ///< 27x5d int vector (stencil for local forces)

typedef std::vector<Vec3s, Eigen::aligned_allocator<Vec3s> >
    Vec3sArray;  ///< an array of 3d float vectors
typedef std::vector<Vec4s, Eigen::aligned_allocator<Vec4s> >
    Vec4sArray;  ///< an array of 4d float vectors
typedef Eigen::Matrix<short, Eigen::Dynamic, Eigen::Dynamic>
    MatXs;  ///< arbitrary dimension scalar matrix

typedef Eigen::Matrix<float, 2, 1> Vec2f;  ///< 2d scalar vector
typedef Eigen::Matrix<float, 3, 1> Vec3f;  ///< 3d scalar vector
typedef Eigen::Matrix<float, 4, 1> Vec4f;  ///< 4d scalar vector
typedef Eigen::Matrix<float, 5, 1> Vec5f;  ///< 4d scalar vector
typedef Eigen::Matrix<float, 6, 1> Vec6f;  ///< 4d scalar vector
typedef Eigen::Matrix<float, 7, 1> Vec7f;  ///< 4d scalar vector
typedef Eigen::Matrix<float, 8, 1> Vec8f;  ///< 4d scalar vector
typedef Eigen::Matrix<float, 9, 1> Vec9f;  ///< 4d scalar vector
typedef Eigen::Matrix<float, 11, 1>
    Vec11f;  ///< 11d scalar vector (stencil for local forces)
typedef Eigen::Matrix<float, 27, 2>
    Mat27x2f;  ///< 27x2d scalar vector (stencil for local forces)
typedef Eigen::Matrix<float, 27, 3>
    Mat27x3f;  ///< 27x3d scalar vector (stencil for local forces)
typedef Eigen::Matrix<float, 27, 4>
    Mat27x4f;  ///< 27x4d scalar vector (stencil for local forces)
typedef Eigen::Matrix<float, 27, 5>
    Mat27x5f;  ///< 27x5d scalar vector (stencil for local forces)
typedef Eigen::Matrix<float, Eigen::Dynamic, 1>
    VecXf;  ///< arbitrary dimension scalar vector

typedef std::vector<Vec3f, Eigen::aligned_allocator<Vec3f> >
    Vec3fArray;  ///< an array of 3d float vectors
typedef std::vector<Vec4f, Eigen::aligned_allocator<Vec4f> >
    Vec4fArray;  ///< an array of 4d float vectors

typedef Eigen::Matrix<double, 2, 1> Vec2d;
typedef Eigen::Matrix<double, 3, 1> Vec3d;
typedef Eigen::Matrix<double, 4, 1> Vec4d;
typedef std::vector<Vec2d, Eigen::aligned_allocator<Vec2d> >
    Vec2dArray;  ///< an array of 2d double vectors
// typedef std::vector<Vec3d, Eigen::aligned_allocator<Vec3d> > Vec3dArray; ///<
// an array of 3d double vectors
typedef std::vector<Vec3d, Eigen::aligned_allocator<Vec3d> >
    Vec3dArray;  ///< an array of 3d double vectors
typedef std::vector<Vec4d, Eigen::aligned_allocator<Vec4d> >
    Vec4dArray;  ///< an array of 3d double vectors
typedef Eigen::Matrix<double, Eigen::Dynamic, 1>
    VecXd;  ///< arbitrary dimension scalar vector

typedef Eigen::Matrix<Scalar, 2, 2> Mat2x;  ///< 2x2 scalar matrix
typedef Eigen::Matrix<Scalar, 3, 3> Mat3x;  ///< 3x3 scalar matrix
typedef Eigen::Matrix<Scalar, 4, 4> Mat4x;  ///< 4x4 scalar matrix
typedef Eigen::Matrix<Scalar, 6, 6> Mat6x;  ///< 4x4 scalar matrix
typedef Eigen::Matrix<Scalar, 11, 11>
    Mat11x;  ///< 11x11 scalar matrix (stencil for local forces)
typedef std::vector<Mat11x, Eigen::aligned_allocator<Mat11x> >
    Mat11xArray;  ///< an array of 11d scalar matrices
typedef std::pair<Mat11x, Mat11x> Mat11xPair;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    MatXx;  ///< arbitrary dimension scalar matrix

typedef Eigen::Matrix<float, 2, 2> Mat2f;  ///< 2x2 scalar matrix
typedef Eigen::Matrix<float, 3, 3> Mat3f;  ///< 3x3 scalar matrix
typedef Eigen::Matrix<float, 4, 4> Mat4f;  ///< 4x4 scalar matrix
typedef Eigen::Matrix<float, 6, 6> Mat6f;  ///< 4x4 scalar matrix
typedef Eigen::Matrix<float, 11, 11>
    Mat11f;  ///< 11x11 scalar matrix (stencil for local forces)
typedef std::vector<Mat11f, Eigen::aligned_allocator<Mat11f> >
    Mat11fArray;  ///< an array of 11d scalar matrices
typedef std::pair<Mat11f, Mat11f> Mat11fPair;
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>
    MatXf;  ///< arbitrary dimension scalar matrix

typedef Eigen::Quaternion<Scalar> Quaternion;

typedef Eigen::SparseMatrix<Scalar, Eigen::ColMajor, int> SparseMatx;
typedef Eigen::SparseMatrix<Scalar, Eigen::RowMajor, int> SparseRowMatx;

template <typename ScalarT>
ScalarT EIGEN_STRONG_INLINE SMALL_NUMBER() {
  return std::numeric_limits<ScalarT>::epsilon();
}

template <>
float EIGEN_STRONG_INLINE SMALL_NUMBER<float>() {
  return 1e-6f;
}

template <>
double EIGEN_STRONG_INLINE SMALL_NUMBER<double>() {
  return 1e-12;
}

EIGEN_STRONG_INLINE Scalar square(const Scalar x) { return x * x; }

EIGEN_STRONG_INLINE Scalar cube(const Scalar x) { return x * x * x; }

template <typename ComparableT>
EIGEN_STRONG_INLINE ComparableT clamp(const ComparableT x, const ComparableT l,
                                      const ComparableT u) {
  return (x > u) ? u : ((x > l) ? x : l);
}

template <typename ScalarT>
EIGEN_STRONG_INLINE bool isSmall(ScalarT x) {
  return fabs(x) < SMALL_NUMBER<ScalarT>();
}

template <typename NormableT>
EIGEN_STRONG_INLINE bool isClose(const NormableT& x1, const NormableT& x2) {
  return isSmall((x1 - x2).norm());
}

template <typename NormableT>
EIGEN_STRONG_INLINE bool isApproxUnit(const NormableT& x) {
  return isSmall(x.squaredNorm() - 1);
}

// Use int as index type ; as there is no MFnUintArrayData in maya API, it
// simplifies conversions There is a MfnUint64ArrayData type though, but that
// would be overkill even if we some day manage to have 150k strands inside the
// simulation
typedef std::vector<int> StrandSet;

}  // namespace strandsim

namespace std {
template <typename Derived>
inline void swap(Eigen::DenseBase<Derived>& a, Eigen::DenseBase<Derived>& b) {
  a.swap(b);
}

template <typename Derived>
inline void swap(
    pair<Eigen::DenseBase<Derived>, Eigen::DenseBase<Derived> >& a,
    pair<Eigen::DenseBase<Derived>, Eigen::DenseBase<Derived> >& b) {
  a.first.swap(b.first);
  a.second.swap(b.second);
}
}  // namespace std

#endif /* S_DEFINITIONS_HH_ */
