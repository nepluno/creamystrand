/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef __MATH_UTILITIES_H__
#define __MATH_UTILITIES_H__

#include <Eigen/Core>
#include <algorithm>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "../Core/Definitions.hh"

#ifndef M_PI
const Scalar M_PI = 3.14159265358979323846264338327950288;
#endif

#ifdef WIN32
#undef min
#undef max
#endif

using std::max;
using std::min;
using std::swap;

namespace strandsim {

/**
 * \brief Matrix representation of the cross product operator.
 */
inline Mat3x crossMat(const Vec3x& a) {
  Mat3x M;
  M << 0, -a(2), a(1), a(2), 0, -a(0), -a(1), a(0), 0;

  return M;
}

/**
 * \brief Parallel-transports u along the t0->t1 transformation.
 *
 * \param [in] u vector to be parallel-transported.
 * \param [in] t0 first vector defining transformation. It doesn't have to be
 * normalized. \param [in] t1 second vector defining transformation. It doesn't
 * have to be normalized.
 *
 * \return If t0 and t1 are colinear, then u is unchanged. Otherwise (t0, t1)
 * defines a plane and a rotation around an axis perpendicular to that plane,
 * which sends t0.normalized() onto t1.normalized(). The parallel transport of u
 * is its image by this rotation.
 *
 * \see normalParallelTransport()
 */

inline Vec3x parallelTransport(const Vec3x& u, const Vec3x& t0,
                               const Vec3x& t1) {
  // Compute rotation axis (if any)
  Vec3x b = t0.cross(t1);
  const Scalar bNorm = b.norm();
  if (isSmall(bNorm))  // vectors are nearly collinear
    return u;
  b /= bNorm;

  const Vec3x& n0 = t0.cross(b).normalized();
  const Vec3x& n1 = t1.cross(b).normalized();

  return u.dot(t0.normalized()) * t1.normalized() + u.dot(n0) * n1 +
         u.dot(b) * b;
}

/**
 * \brief Parallel-transports u along the t0->t1 transformation, assuming that u
 * is normal to t0.
 *
 * \param [in] u vector to be parallel-transported.
 * \param [in] t0 first vector defining transformation. It doesn't have to be
 * normalized. \param [in] t1 second vector defining transformation. It doesn't
 * have to be normalized.
 *
 * \return If t0 and t1 are colinear, then u is unchanged. Otherwise (t0, t1)
 * defines a plane and a rotation around an axis perpendicular to that plane,
 * which sends t0.normalized() onto t1.normalized(). The parallel transport of u
 * is its image by this rotation.
 *
 * \note It is assumed (see assertion below) that u.dot(0)==0, which is the case
 * when parallel-transporting frame vectors for instance. Then this is cheaper
 * to call than parallelTransport()
 */
inline Vec3x normalParallelTransport(const Vec3x& u, const Vec3x& t0,
                                     const Vec3x& t1) {
  // This should be called only to transport an orthogonal vector
  assert(isSmall(u.dot(t0)));

  // Compute rotation axis (if any)
  Vec3x b = t0.cross(t1);
  const Scalar bNorm = b.norm();
  if (isSmall(bNorm))  // vectors are nearly collinear
    return u;
  b /= bNorm;

  const Vec3x& n0 = t0.cross(b).normalized();
  const Vec3x& n1 = t1.cross(b).normalized();

  return u.dot(n0) * n1 + u.dot(b) * b;
}

/**
 * \brief Parallel-transports u along the t0->t1 transformation, assuming that u
 * is normal to t0. and all the vectors are unitary
 *
 * \param [in] u vector to be parallel-transported.
 * \param [in] t0 first vector defining transformation. It  have to be
 * normalized. \param [in] t1 second vector defining transformation. It have to
 * be normalized.
 *
 * \return If t0 and t1 are colinear, then u is unchanged. Otherwise (t0, t1)
 * defines a plane and a rotation around an axis perpendicular to that plane,
 * which sends t0 onto t1. The parallel transport of u is its image by this
 * rotation.
 *
 * \note It is assumed that u.dot( t0 ) = 0 and t0.norm() = t1.norm() = u.norm()
 * = 1
 */
inline Vec3x orthonormalParallelTransport(const Vec3x& u, const Vec3x& t0,
                                          const Vec3x& t1) {
  // This should be called only to transport an orthogonal vector
  assert(isSmall(u.dot(t0)));

  Vec3x b = t0.cross(t1);
  const Scalar bNorm = b.norm();
  if (isSmall(bNorm))  // vectors are nearly collinear
    return u;
  b /= bNorm;

  const Vec3x& n0 = t0.cross(b);
  const Vec3x& n1 = t1.cross(b);

  return u.dot(n0) * n1 + u.dot(b) * b;
}

inline bool containsNans(const VecXx& dofs) {
  bool panic = false;

  for (int i = 0; i < dofs.size(); ++i) {
    if (std::isnan(dofs[i])) {
      panic = true;
      break;
    }
  }

  return panic;
}

template <typename T, typename P>
inline void check_isnan(const char* name, const T& val, const P& toprint) {
  if (std::isnan(val)) {
    std::cout << "NAN in " << name << std::endl;
    std::cout << toprint << std::endl;
    exit(EXIT_FAILURE);
  } else if (val > 1e+5 || val < -1e+5) {
    std::cout << "WARNING: LARGE " << name << std::endl;
    std::cout << toprint << std::endl;

    if (std::isinf(val)) {
      std::cout << "INF " << name << std::endl;
      std::cout << toprint << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

template <typename T>
inline void check_isnan(const char* name, const T& val) {
  check_isnan(name, val, val);
}

template <typename T>
inline void check_isnan(const char* name, const std::vector<T>& v) {
  int i = 0;
  for (const T& val : v) {
    std::ostringstream oss;
    oss << name << " [" << i << "]";
    check_isnan(oss.str().c_str(), val);
    ++i;
  }
}

template <typename T>
inline void check_isnan(
    const char* name,
    const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& v) {
  int i = 0;
  for (const auto& vec : v) {
    std::ostringstream oss;
    oss << name << " [" << i << "]";
    check_isnan(oss.str().c_str(),
                sqrt(vec.squaredNorm() / (double)std::max(1, (int)vec.size())),
                vec);
    ++i;
  }
}

template <typename T>
inline void check_isnan(
    const char* name,
    const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vx,
    const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vy,
    const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vz) {
  std::string name0(name);
  name0 += ", x";
  check_isnan(name0.c_str(), vx);

  std::string name1(name);
  name1 += ", y";
  check_isnan(name1.c_str(), vy);

  std::string name2(name);
  name2 += ", z";
  check_isnan(name2.c_str(), vz);
}

template <typename T, int M, int N>
inline void check_isnan(const char* name, const Eigen::Matrix<T, M, N>& v) {
  check_isnan(name, sqrt(v.squaredNorm() / (double)std::max(1, (int)v.size())),
              v);
}

template <typename T>
inline void check_isnan(const char* name, int N, const T& func) {
  for (int i = 0; i < N; ++i) {
    check_isnan(name, func(i));
  }
}

template <class S, class T>
inline S lerp(const S& value0, const S& value1, T f) {
  return (1 - f) * value0 + f * value1;
}

template <class S, class T>
inline S bilerp(const S& v00, const S& v10, const S& v01, const S& v11, T fx,
                T fy) {
  return lerp(lerp(v00, v10, fx), lerp(v01, v11, fx), fy);
}

template <class S, class T>
inline S trilerp(const S& v000, const S& v100, const S& v010, const S& v110,
                 const S& v001, const S& v101, const S& v011, const S& v111,
                 T fx, T fy, T fz) {
  return lerp(bilerp(v000, v100, v010, v110, fx, fy),
              bilerp(v001, v101, v011, v111, fx, fy), fz);
}

template <class T>
void zero(std::vector<T>& v) {
  for (int i = (int)v.size() - 1; i >= 0; --i) v[i] = 0;
}

template <class T>
T abs_max(const std::vector<T>& v) {
  T m = 0;
  for (int i = (int)v.size() - 1; i >= 0; --i) {
    if (std::fabs(v[i]) > m) m = std::fabs(v[i]);
  }
  return m;
}

template <class T>
bool contains(const std::vector<T>& a, T e) {
  for (unsigned int i = 0; i < a.size(); ++i)
    if (a[i] == e) return true;
  return false;
}

template <class T>
void add_unique(std::vector<T>& a, T e) {
  for (unsigned int i = 0; i < a.size(); ++i)
    if (a[i] == e) return;
  a.push_back(e);
}

template <class T>
void insert(std::vector<T>& a, unsigned int index, T e) {
  a.push_back(a.back());
  for (unsigned int i = (unsigned int)a.size() - 1; i > index; --i)
    a[i] = a[i - 1];
  a[index] = e;
}
template <class T>
void erase(std::vector<T>& a, unsigned int index) {
  for (unsigned int i = index; i < a.size() - 1; ++i) a[i] = a[i + 1];
  a.pop_back();
}

template <class T>
void erase_swap(std::vector<T>& a, unsigned int index) {
  for (unsigned int i = index; i < a.size() - 1; ++i) swap(a[i], a[i + 1]);
  a.pop_back();
}

template <class T>
void erase_unordered(std::vector<T>& a, unsigned int index) {
  a[index] = a.back();
  a.pop_back();
}

template <class T>
void erase_unordered_swap(std::vector<T>& a, unsigned int index) {
  swap(a[index], a.back());
  a.pop_back();
}

template <class T>
void find_and_erase_unordered(std::vector<T>& a, const T& doomed_element) {
  for (unsigned int i = 0; i < a.size(); ++i)
    if (a[i] == doomed_element) {
      erase_unordered(a, i);
      return;
    }
}

template <class T>
void replace_once(std::vector<T>& a, const T& old_element,
                  const T& new_element) {
  for (unsigned int i = 0; i < a.size(); ++i)
    if (a[i] == old_element) {
      a[i] = new_element;
      return;
    }
}
template <class T>
inline Eigen::Matrix<T, 3, 1> grad_trilerp(const T& v000, const T& v100,
                                           const T& v010, const T& v110,
                                           const T& v001, const T& v101,
                                           const T& v011, const T& v111, T fx,
                                           T fy, T fz) {
  return Eigen::Matrix<T, 3, 1>(-(fy - 1.) * (fz - 1.), -(fx - 1.) * (fz - 1.),
                                -(fx - 1.) * (fy - 1.)) *
             v000 +
         Eigen::Matrix<T, 3, 1>((fy - 1.) * (fz - 1.), fx * (fz - 1.),
                                fx * (fy - 1.)) *
             v100 +
         Eigen::Matrix<T, 3, 1>(fy * (fz - 1.), (fx - 1.) * (fz - 1.),
                                fy * (fx - 1.)) *
             v010 +
         Eigen::Matrix<T, 3, 1>(-fy * (fz - 1.), -fx * (fz - 1.), -fx * fy) *
             v110 +
         Eigen::Matrix<T, 3, 1>(fz * (fy - 1.), fz * (fx - 1.),
                                (fx - 1.) * (fy - 1.)) *
             v001 +
         Eigen::Matrix<T, 3, 1>(-fz * (fy - 1.), -fx * fz, -fx * (fy - 1.)) *
             v101 +
         Eigen::Matrix<T, 3, 1>(-fy * fz, -fz * (fx - 1.), -fy * (fx - 1.)) *
             v011 +
         Eigen::Matrix<T, 3, 1>(fy * fz, fx * fz, fx * fy) * v111;
}

template <class T>
inline void get_barycentric(T x, int& i, T& f, int i_low, int i_high) {
  T s = std::floor(x);
  i = (int)s;
  if (i < i_low) {
    i = i_low;
    f = 0;
  } else if (i > i_high - 2) {
    i = i_high - 2;
    f = 1;
  } else
    f = (T)(x - s);
}

template <typename S>
inline void make_gibbs_simplex(Eigen::Matrix<S, Eigen::Dynamic, 1>& v) {
  S vs = v.sum();
  if (vs < (S)1e-20) {
    v.setZero();
    v(0) = (S)1.0;
  } else {
    v /= vs;
  }
}

template <typename S, int K>
inline void swap(Eigen::Matrix<S, Eigen::Dynamic, 1>& v, int i, int j) {
  Eigen::Matrix<S, K, 1> c = v.template segment<K>(i * K);
  v.template segment<K>(i * K) = v.template segment<K>(j * K);
  v.template segment<K>(j * K) = c;
}

template <typename S>
inline void swap(Eigen::Matrix<S, Eigen::Dynamic, 1>& v, int i, int j, int K) {
  Eigen::Matrix<S, Eigen::Dynamic, 1> c = v.segment(i * K, K);
  v.segment(i * K, K) = v.segment(j * K, K);
  v.segment(j * K, K) = c;
}

template <typename S, int K>
inline void swap(Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>& v, int i,
                 int j) {
  Eigen::Matrix<S, K, K> c = v.template block<K, K>(i * K, 0);
  v.template block<K, K>(i * K, 0) = v.template block<K, K>(j * K, 0);
  v.template block<K, K>(j * K, 0) = c;
}

inline void fisherYates(int n, std::vector<int>& indices) {
  indices.resize(n);
  for (int i = 0; i < n; ++i) indices[i] = i;
  for (int i = n - 1; i >= 1; --i) {
    const int j = rand() % (i + 1);
    std::swap(indices[i], indices[j]);
  }
}

bool approxSymmetric(const MatXx& A, const Scalar& eps);

Scalar ScalarRand(const Scalar min, const Scalar max);

inline Scalar perimeter(const Scalar& ra, const Scalar& rb) {
  // The second-Ramanujan approximation to ellipse perimeter
  const Scalar h = (ra - rb) * (ra - rb) / ((ra + rb) * (ra + rb));
  return M_PI * (ra + rb) * (1. + 3. * h / (10. + sqrt(4. - 3. * h)));
}

inline int mod_floor(int a, int n) { return ((a % n) + n) % n; }

inline Scalar sgn(const Scalar& x) {
  return x == 0.0 ? 0.0 : (x > 0.0 ? 1.0 : -1.0);
}

inline Scalar mean_curvature(const Vec9x& h, const Scalar& dx) {
  const Scalar t2 = 1.0 / (dx * dx);
  const Scalar t3 = h[1] - h[7];
  const Scalar t7 = h[4] * 2.0;
  const Scalar t4 = h[3] + h[5] - t7;
  const Scalar t5 = 1.0 / (dx * dx * dx * dx);
  const Scalar t6 = h[3] - h[5];
  const Scalar t8 = dx * dx;
  const Scalar t0 =
      pow(t8, 3.0 / 2.0) *
      (t2 * t4 + t2 * (h[1] - h[4] * 2.0 + h[7]) +
       (t3 * t3) * t4 * t5 * (1.0 / 4.0) +
       t5 * (t6 * t6) * (h[1] + h[7] - t7) * (1.0 / 4.0) -
       t3 * t5 * t6 * (h[0] - h[2] - h[6] + h[8]) * (1.0 / 8.0)) *
      1.0 /
      pow(t8 * 4.0 - h[1] * h[7] * 2.0 - h[3] * h[5] * 2.0 + h[1] * h[1] +
              h[3] * h[3] + h[5] * h[5] + h[7] * h[7],
          3.0 / 2.0) *
      8.0;
  return t0;
}

inline Scalar wedge_wet_abscissa(const Scalar& A0, const Scalar& A1,
                                 const Scalar& d0, const Scalar& d1,
                                 const Scalar& theta) {
  if (fabs(d0 - d1) < 1e-12) {
    if (fabs(A0 - A1) < 1e-12) {
      if (d0 < (1. + 0.5 * theta) * sqrt(A0)) {
        return 1.0;
      } else {
        return 0.0;
      }
    } else {
      return clamp(
          (A0 - (d0 * d0) * 1.0 / pow(theta * (1.0 / 2.0) + 1.0, 2.0)) /
              (A0 - A1),
          0.0, 1.0);
    }
  } else {
    const Scalar t2 = theta * theta;
    const Scalar t3 = A0 * A0;
    const Scalar t4 = A1 * A1;
    const Scalar t5 = d0 * d0;
    const Scalar t6 = d1 * d1;
    const Scalar t7 = A0 * t6 * 1.6E1;
    const Scalar t8 = A1 * t5 * 1.6E1;
    const Scalar t9 = t3 * theta * 4.0;
    const Scalar t10 = t4 * theta * 4.0;
    const Scalar t11 = t3 * 4.0;
    const Scalar t12 = t4 * 4.0;
    const Scalar t13 = t2 * t3;
    const Scalar t14 = t2 * t4;
    const Scalar t15 = t7 + t8 + t9 + t10 + t11 + t12 + t13 + t14 -
                       A0 * A1 * 8.0 - A0 * A1 * t2 * 2.0 -
                       A0 * A1 * theta * 8.0 - A0 * d0 * d1 * 1.6E1 -
                       A1 * d0 * d1 * 1.6E1;
    if (t15 < 0.0) {
      return 0.0;
    }
    const Scalar t16 = sqrt(t15);
    return clamp(
        1.0 / pow(d0 - d1, 2.0) *
            (A0 * -4.0 + A1 * 4.0 + t5 * 8.0 + t16 * 2.0 - A0 * t2 + A1 * t2 -
             A0 * theta * 4.0 + A1 * theta * 4.0 - d0 * d1 * 8.0 + t16 * theta),
        0.0, 1.0);
  }
}

// Given two signed distance values (line endpoints), determine what fraction of
// a connecting segment is "inside"
inline Scalar fraction_inside(const Scalar& phi_left, const Scalar& phi_right) {
  if (phi_left < 0 && phi_right < 0) return 1;
  if (phi_left < 0 && phi_right >= 0) return phi_left / (phi_left - phi_right);
  if (phi_left >= 0 && phi_right < 0)
    return phi_right / (phi_right - phi_left);
  else
    return 0;
}

inline void cycle_array(Scalar* arr, int size) {
  Scalar t = arr[0];
  for (int i = 0; i < size - 1; ++i) arr[i] = arr[i + 1];
  arr[size - 1] = t;
}

inline void d2WsdF2B(const Scalar J, const Mat3x& F, const Mat3x& b_bar,
                     const Scalar& mu, const Mat3x& B, Mat3x& ret) {
  Mat3x invF = F.inverse();
  Mat3x devb_bar = b_bar - 1.0 / 3.0 * b_bar.trace() * Mat3x::Identity();
  Mat3x invFT = invF.transpose();
  Mat3x M0 = devb_bar * invFT;
  Scalar J_fac = pow(J, -2.0 / 3.0);

  ret = mu * (J_fac * B - 2.0 / 3.0 * (invF * B).trace() * M0 -
              1.0 / 3.0 * invFT *
                  (J_fac * (F.transpose() * B).trace() * Mat3x::Identity() -
                   b_bar.trace() * B.transpose() * invFT));
}

inline Scalar fraction_inside(const Scalar& phi_bl, const Scalar& phi_br,
                              const Scalar& phi_tl, const Scalar& phi_tr) {
  int inside_count = (phi_bl < 0 ? 1 : 0) + (phi_tl < 0 ? 1 : 0) +
                     (phi_br < 0 ? 1 : 0) + (phi_tr < 0 ? 1 : 0);
  Scalar list[] = {phi_bl, phi_br, phi_tr, phi_tl};

  if (inside_count == 4)
    return 1.;
  else if (inside_count == 3) {
    // rotate until the positive value is in the first position
    while (list[0] < 0) {
      cycle_array(list, 4);
    }

    // Work out the area of the exterior triangle
    Scalar side0 = 1. - fraction_inside(list[0], list[3]);
    Scalar side1 = 1. - fraction_inside(list[0], list[1]);
    return 1. - 0.5 * side0 * side1;
  } else if (inside_count == 2) {
    // rotate until a negative value is in the first position, and the next
    // negative is in either slot 1 or 2.
    while (list[0] >= 0 || !(list[1] < 0 || list[2] < 0)) {
      cycle_array(list, 4);
    }

    if (list[1] < 0) {  // the matching signs are adjacent
      Scalar side_left = fraction_inside(list[0], list[3]);
      Scalar side_right = fraction_inside(list[1], list[2]);
      return 0.5 * (side_left + side_right);
    } else {  // matching signs are diagonally opposite
      // determine the centre point's sign to disambiguate this case
      Scalar middle_point = 0.25 * (list[0] + list[1] + list[2] + list[3]);
      if (middle_point < 0) {
        Scalar area = 0.;

        // first triangle (top left)
        Scalar side1 = 1. - fraction_inside(list[0], list[3]);
        Scalar side3 = 1. - fraction_inside(list[2], list[3]);

        area += 0.5 * side1 * side3;

        // second triangle (top right)
        Scalar side2 = 1 - fraction_inside(list[2], list[1]);
        Scalar side0 = 1 - fraction_inside(list[0], list[1]);
        area += 0.5 * side0 * side2;

        return 1. - area;
      } else {
        Scalar area = 0.;

        // first triangle (bottom left)
        Scalar side0 = fraction_inside(list[0], list[1]);
        Scalar side1 = fraction_inside(list[0], list[3]);
        area += 0.5 * side0 * side1;

        // second triangle (top right)
        Scalar side2 = fraction_inside(list[2], list[1]);
        Scalar side3 = fraction_inside(list[2], list[3]);
        area += 0.5 * side2 * side3;
        return area;
      }
    }
  } else if (inside_count == 1) {
    // rotate until the negative value is in the first position
    while (list[0] >= 0) {
      cycle_array(list, 4);
    }

    // Work out the area of the interior triangle, and subtract from 1.
    Scalar side0 = fraction_inside(list[0], list[3]);
    Scalar side1 = fraction_inside(list[0], list[1]);
    return 0.5 * side0 * side1;
  } else
    return 0;
}

template <typename S, unsigned N, unsigned M>
inline Eigen::Matrix<S, N, M> frac(const Eigen::Matrix<S, N, M>& x) {
  Eigen::Matrix<S, N, M> ret(x);
  for (unsigned i = 0; i < N; ++i)
    for (unsigned j = 0; j < M; ++j) ret(i, j) -= floor(ret(i, j));
  return ret;
}

inline Scalar sqr(const Scalar& s) { return s * s; }

inline Vec3x cross_x(const Vec3x& a, const Scalar& b) {
  return Vec3x(0.0, a(2) * b, -a(1) * b);
}

inline Vec3x cross_y(const Vec3x& a, const Scalar& b) {
  return Vec3x(-a(2) * b, 0.0, a(0) * b);
}

inline Vec3x cross_z(const Vec3x& a, const Scalar& b) {
  return Vec3x(a(1) * b, -a(0) * b, 0.0);
}

inline Scalar cross_x_row(const Vec3x& a, const Vec3x& b) {
  return -a(2) * b(1) + a(1) * b(2);
}

inline Scalar cross_y_row(const Vec3x& a, const Vec3x& b) {
  return a(2) * b(0) - a(0) * b(2);
}

inline Scalar cross_z_row(const Vec3x& a, const Vec3x& b) {
  return -a(1) * b(0) + a(0) * b(1);
}

template <typename T>
inline T smooth_kernel(const T& r2, const T& h) {
  return std::max((T)pow((T)1.0 - r2 / (h * h), (T)3.0), (T)0.0);
}

template <typename T>
inline T grad_smooth_kernel_dot_dist(const T& r2, const T& h) {
  T r2h2 = std::min(1.0, r2 / (h * h));
  return ((T)1.0 - r2h2) * ((T)1.0 - r2h2) * r2h2;
}

template <typename T>
inline T poly6_kernel(const T& r2, const T& h) {
  return std::max((T)pow(h * h - r2, (T)3.0), (T)0.0) / pow(h, 9.0) *
         1.5666814711;  // unit: cm^{-3}
}

inline Scalar quad_kernel(const Scalar& x) {
  const Scalar dx = fabs(x);
  if (dx < 0.5)
    return 0.75 - dx * dx;
  else if (dx < 1.5)
    return 0.5 * (1.5 - dx) * (1.5 - dx);
  else
    return 0.0;
}

inline Scalar cubic_kernel(const Scalar& x) {
  const Scalar dx = fabs(x);
  if (dx < 1.0)
    return 0.5 * dx * dx * dx - dx * dx + 2.0 / 3.0;
  else if (dx < 2.0)
    return (2.0 - dx) * (2.0 - dx) * (2.0 - dx) / 6.0;
  else
    return 0.0;
}

inline Scalar linear_kernel(const Scalar& x) {
  const Scalar dx = fabs(x);
  if (dx < 1.0)
    return 1.0 - dx;
  else
    return 0.0;
}

inline Scalar grad_quad_kernel(const Scalar& x) {
  const Scalar dx = fabs(x);
  if (dx < 0.5)
    return -2.0 * dx * sgn(x);
  else if (dx < 1.5)
    return sgn(x) * (dx - 1.5);
  else
    return 0.0;
}

inline Scalar grad_cubic_kernel(const Scalar& x) {
  const Scalar dx = fabs(x);
  if (dx < 1.0)
    return (x * sgn(x) * (3.0 * x - 4.0 * sgn(x))) / 2.0;
  else if (dx < 2.0)
    return -(sgn(x) * (dx - 2.0) * (dx - 2.0)) / 2.0;
  else
    return 0.0;
}

inline Scalar grad_linear_kernel(const Scalar& x) {
  const Scalar dx = fabs(x);
  if (dx < 1.0)
    return -sgn(x);
  else
    return 0.0;
}

template <int order>
inline Scalar N_kernel(const Vec3x& x) {
  if (order == 2)
    return quad_kernel(x(0)) * quad_kernel(x(1)) * quad_kernel(x(2));
  else if (order == 3)
    return cubic_kernel(x(0)) * cubic_kernel(x(1)) * cubic_kernel(x(2));
  else
    return linear_kernel(x(0)) * linear_kernel(x(1)) * linear_kernel(x(2));
}

template <int order>
inline Scalar grad_N_kernel(const Vec3x& x, const Scalar& h, Vec3x& g) {
  if (order == 1) {
    const Scalar v0 = linear_kernel(x(0));
    const Scalar v1 = linear_kernel(x(1));
    const Scalar v2 = linear_kernel(x(2));
    const Scalar g0 = grad_linear_kernel(x(0));
    const Scalar g1 = grad_linear_kernel(x(1));
    const Scalar g2 = grad_linear_kernel(x(2));

    g(0) = g0 * v1 * v2;
    g(1) = v0 * g1 * v2;
    g(2) = v0 * v1 * g2;
    g *= 1.0 / h;
    return v0 * v1 * v2;
  } else if (order == 2) {
    const Scalar inv_h = 1. / h;
    const Scalar dx0 = fabs(x(0));
    const Scalar dx1 = fabs(x(1));
    const Scalar dx2 = fabs(x(2));
    const Scalar sx0 = sgn(x(0));
    const Scalar sx1 = sgn(x(1));
    const Scalar sx2 = sgn(x(2));

    const Scalar v0 =
        dx0 < 0.5 ? (0.75 - dx0 * dx0)
                  : (dx0 < 1.5 ? (0.5 * (1.5 - dx0) * (1.5 - dx0)) : 0.0);
    const Scalar v1 =
        dx1 < 0.5 ? (0.75 - dx1 * dx1)
                  : (dx1 < 1.5 ? (0.5 * (1.5 - dx1) * (1.5 - dx1)) : 0.0);
    const Scalar v2 =
        dx2 < 0.5 ? (0.75 - dx2 * dx2)
                  : (dx2 < 1.5 ? (0.5 * (1.5 - dx2) * (1.5 - dx2)) : 0.0);

    const Scalar g0 =
        dx0 < 0.5 ? (-2.0 * dx0 * sx0) : (dx0 < 1.5 ? sx0 * (dx0 - 1.5) : 0.0);
    const Scalar g1 =
        dx1 < 0.5 ? (-2.0 * dx1 * sx1) : (dx1 < 1.5 ? sx1 * (dx1 - 1.5) : 0.0);
    const Scalar g2 =
        dx2 < 0.5 ? (-2.0 * dx2 * sx2) : (dx2 < 1.5 ? sx2 * (dx2 - 1.5) : 0.0);

    g(0) = g0 * v1 * v2 * inv_h;
    g(1) = v0 * g1 * v2 * inv_h;
    g(2) = v0 * v1 * g2 * inv_h;
    return v0 * v1 * v2;
  } else {
    const Scalar v0 = cubic_kernel(x(0));
    const Scalar v1 = cubic_kernel(x(1));
    const Scalar v2 = cubic_kernel(x(2));
    const Scalar g0 = grad_cubic_kernel(x(0));
    const Scalar g1 = grad_cubic_kernel(x(1));
    const Scalar g2 = grad_cubic_kernel(x(2));

    g(0) = g0 * v1 * v2;
    g(1) = v0 * g1 * v2;
    g(2) = v0 * v1 * g2;
    g *= 1.0 / h;
    return v0 * v1 * v2;
  }
}

inline Scalar grad_N_kernel_x(const Vec3x& x, const Scalar& h,
                              const int order) {
  std::function<Scalar(const Scalar&)> val =
      (order == 1) ? linear_kernel
                   : ((order == 3) ? cubic_kernel : quad_kernel);
  std::function<Scalar(const Scalar&)> grad =
      (order == 1) ? grad_linear_kernel
                   : ((order == 3) ? grad_cubic_kernel : grad_quad_kernel);

  return grad(x(0)) * val(x(1)) * val(x(2)) / h;
}

inline Scalar grad_N_kernel_y(const Vec3x& x, const Scalar& h,
                              const int order) {
  std::function<Scalar(const Scalar&)> val =
      (order == 1) ? linear_kernel
                   : ((order == 3) ? cubic_kernel : quad_kernel);
  std::function<Scalar(const Scalar&)> grad =
      (order == 1) ? grad_linear_kernel
                   : ((order == 3) ? grad_cubic_kernel : grad_quad_kernel);

  return val(x(0)) * grad(x(1)) * val(x(2)) / h;
}

inline Scalar grad_N_kernel_z(const Vec3x& x, const Scalar& h,
                              const int order) {
  std::function<Scalar(const Scalar&)> val =
      (order == 1) ? linear_kernel
                   : ((order == 3) ? cubic_kernel : quad_kernel);
  std::function<Scalar(const Scalar&)> grad =
      (order == 1) ? grad_linear_kernel
                   : ((order == 3) ? grad_cubic_kernel : grad_quad_kernel);

  return val(x(0)) * val(x(1)) * grad(x(2)) / h;
}

template <typename S, unsigned N>
inline void QRDecompose(const Eigen::Matrix<S, N, N>& A,
                        Eigen::Matrix<S, N, N>& Q, Eigen::Matrix<S, N, N>& R) {
  Eigen::HouseholderQR<Eigen::Matrix<S, N, N> > qr(A);
  R = qr.matrixQR().template triangularView<Eigen::Upper>();
  VecXx diag = VecXx(R.diagonal().array());
  for (int i = 0; i < diag.size(); ++i) diag(i) = sgn(diag(i));
  MatXx s = diag.asDiagonal();
  //        std::cout<<"R: "<<R<<std::endl;
  //        std::cout<<"s:"<<s<<std::endl;
  Q = qr.householderQ() * s;
  R = s * R;
}

inline Scalar defaultRadiusMultiplier() {
  return 0.6203504909;  // pow(0.75 / pi, 1. / 3.)
}

template <typename S, typename T>
inline S lerp_weno(const S value[], T f) {
  S p1 = value[0] + (value[1] - value[0]) * (f + 2.0) +
         (value[2] - 2.0 * value[1] + value[0]) * (f + 2.0) * (f + 1.0) * 0.5 +
         (value[3] - 3.0 * value[2] + 3.0 * value[1] - value[0]) * (f + 2.0) *
             (f + 1.0) * f / 6.0;
  S p2 = value[1] + (value[2] - value[1]) * (f + 1.0) +
         (value[3] - 2.0 * value[2] + value[1]) * (f + 1.0) * f * 0.5 +
         (value[4] - 3.0 * value[3] + 3.0 * value[2] - value[1]) * (f + 1.0) *
             f * (f - 1.0) / 6.0;
  S p3 = value[2] + (value[3] - value[2]) * f +
         (value[4] - 2.0 * value[3] + value[2]) * f * (f - 1.0) * 0.5 +
         (value[5] - 3.0 * value[4] + 3.0 * value[3] - value[2]) * f *
             (f - 1.0) * (f - 2.0) / 6.0;

  T C1 = (2 - f) * (3 - f) / 20.0;
  T C2 = (3 - f) * (f + 2) / 10.0;
  T C3 = (f + 2) * (f + 1) / 20.0;

  T IS1 = (814.0 * value[3] * value[3] + 4326 * value[2] * value[2] +
           2976 * value[1] * value[1] + 244 * value[0] * value[0] -
           3579 * value[2] * value[3] - 6927 * value[2] * value[1] +
           1854 * value[2] * value[0] + 2634 * value[3] * value[1] -
           683 * value[3] * value[0] - 1659 * value[1] * value[0]) /
          180.0;
  T IS2 = (1986 * value[3] * value[3] + 1986 * value[2] * value[2] +
           244 * value[1] * value[1] + 244 * value[4] * value[4] +
           1074 * value[2] * value[4] - 3777 * value[2] * value[3] -
           1269 * value[2] * value[1] + 1074 * value[3] * value[1] -
           1269 * value[4] * value[3] - 293 * value[4] * value[1]) /
          180.0;
  T IS3 = (814 * value[2] * value[2] + 4326 * value[3] * value[3] +
           2976 * value[4] * value[4] + 244 * value[5] * value[5] -
           683 * value[2] * value[5] + 2634 * value[2] * value[4] -
           3579 * value[2] * value[3] - 6927 * value[3] * value[4] +
           1854 * value[3] * value[5] - 1659 * value[4] * value[5]) /
          180.0;

  const T epsilon = 1e-6;
  T alpha1 = C1 / ((IS1 + epsilon) * (IS1 + epsilon));
  T alpha2 = C2 / ((IS2 + epsilon) * (IS2 + epsilon));
  T alpha3 = C3 / ((IS3 + epsilon) * (IS3 + epsilon));

  T sumalpha = alpha1 + alpha2 + alpha3;
  T w1 = alpha1 / sumalpha;
  T w2 = alpha2 / sumalpha;
  T w3 = alpha3 / sumalpha;

  return p1 * w1 + p2 * w2 + p3 * w3;
}

inline Scalar cosine_ease_function(const Scalar& t, const Scalar& t0,
                                   const Scalar& t1, const Scalar& ta,
                                   const Scalar& tb, const Scalar& amp,
                                   const Scalar& freq) {
  const Scalar Ta = ta - t0;
  const Scalar Tt = t - t0;
  const Scalar Tb = t1 - tb;
  const Scalar Te = t1 - t0;
  const Scalar w = 2.0 * M_PI * freq;

  if (t < t0 || t > t1)
    return 0.0;
  else if (t < ta) {
    const Scalar t2 = Ta * w;
    const Scalar t3 = cos(t2);
    const Scalar t4 = sin(t2);
    return 1.0 / (Ta * Ta * Ta) * (Tt * Tt) *
               (amp * t3 * 2.0 + amp * Ta * t4 * w) * -3.0 +
           amp * 1.0 / (Ta * Ta) * Tt * (t3 * 3.0 + Ta * t4 * w) * 2.0;
  } else if (t > tb) {
    const Scalar t2 = Tb - Te;
    const Scalar t3 = t2 * w;
    const Scalar t4 = Te - Tt;
    const Scalar t5 = 1.0 / (Tb * Tb * Tb);
    const Scalar t6 = cos(t3);
    const Scalar t7 = sin(t3);
    return -amp * t5 * (Te * 2.0 - Tt * 2.0) *
               (Tb * t6 * 3.0 - Te * t6 * 2.0 + t6 * Tt * 2.0 +
                (Tb * Tb) * t7 * w + Tb * t7 * w * Tt - Tb * Te * t7 * w) +
           amp * (t4 * t4) * t5 * (t6 * 2.0 + Tb * t7 * w);
  } else {
    return -amp * w * sin(w * Tt);
  }
}

inline Scalar weno_ease_function(const Scalar& t, const Scalar& dt,
                                 const Scalar& t0, const Scalar& t1,
                                 const Scalar& base_dt, const Scalar& cur_pos,
                                 const std::vector<Scalar>& bases) {
  if (t < t0 || t > t1) return 0.0;

  Scalar spos = (t - t0) / base_dt;
  int ipos = (int)spos;
  Scalar fpos = spos - (Scalar)ipos;

  int nb = (int)bases.size();

  Scalar extracted_bases[] = {
      bases[clamp(ipos - 2, 0, nb - 1)], bases[clamp(ipos - 1, 0, nb - 1)],
      bases[clamp(ipos - 0, 0, nb - 1)], bases[clamp(ipos + 1, 0, nb - 1)],
      bases[clamp(ipos + 2, 0, nb - 1)], bases[clamp(ipos + 3, 0, nb - 1)]};

  Scalar target_pos = lerp_weno(extracted_bases, fpos);
  return (target_pos - cur_pos) / dt;
}

inline Scalar cubic_ease_function(const Scalar& t, const Scalar& t0,
                                  const Scalar& t1, const Scalar& ta,
                                  const Scalar& tb, const Scalar& L) {
  Scalar yh = (L * 2.0) / (t1 - t0 + tb - ta);
  if (t < t0 || t > t1)
    return 0.0;
  else {
    if (t < ta)
      return (yh * (t0 - t) * (t0 - t) * (t0 - 3.0 * ta + 2.0 * t)) /
             ((t0 - ta) * (t0 - ta) * (t0 - ta));
    else if (t > tb)
      return (yh * (t1 - t) * (t1 - t) * (t1 - 3.0 * tb + 2.0 * t)) /
             ((t1 - tb) * (t1 - tb) * (t1 - tb));
    else
      return yh;
  }
}

inline Scalar inverse_D_coeff(const Scalar& h, const int order) {
  switch (order) {
    case 1:
      return 6.0 / (h * h);
    case 3:
      return 3.0 / (h * h);
    default:
      return 4.0 / (h * h);
  }
}

inline Scalar D_coeff(const Scalar& h, const int order) {
  switch (order) {
    case 1:
      return (h * h) / 6.0;
    case 3:
      return (h * h) / 3.0;
    default:
      return (h * h) / 4.0;
  }
}

// [Botsch et al. 2010] (3.9)
inline void grad_triangle(const Vec3x& x0, const Vec3x& x1, const Vec3x& x2,
                          Mat3x& grad_coeff) {
  Vec3x n = (x2 - x0).cross(x1 - x0);
  const Scalar len = n.norm();
  if (len < 1e-63) {
    grad_coeff.setZero();
    return;
  }

  n /= len;

  Mat3x mrot = Eigen::AngleAxis<Scalar>(-M_PI / 2.0, n).toRotationMatrix();
  Vec3x v0 = mrot * (x0 - x2) / len;
  Vec3x v1 = mrot * (x1 - x0) / len;

  grad_coeff.block<3, 1>(0, 0) = -(v0 + v1);
  grad_coeff.block<3, 1>(0, 1) = v0;
  grad_coeff.block<3, 1>(0, 2) = v1;
}

inline void get_div_triangle(const Scalar& va0, const Scalar& va1,
                             const Scalar& va2, const Scalar& thickness0,
                             const Scalar& thickness1, const Scalar& thickness2,
                             const Vec3x& x0, const Vec3x& x1, const Vec3x& x2,
                             Vec3x& div0, Vec3x& div1, Vec3x& div2) {
  Vec3x n = (x2 - x0).cross(x1 - x0);
  const Scalar len = n.norm();
  if (len < 1e-63) {
    return;
  }

  n /= len;

  Mat3x mrot = Eigen::AngleAxis<Scalar>(-M_PI / 2.0, n).toRotationMatrix();

  div0 = -mrot * (x2 - x1) * 0.5 * thickness0 / va0;
  div1 = -mrot * (x0 - x2) * 0.5 * thickness1 / va1;
  div2 = -mrot * (x1 - x0) * 0.5 * thickness2 / va2;
}

inline Scalar get_rotated_pore_coeff_x(const Vec3x& a, const Scalar& sf) {
  const Scalar t2 = a(0) * a(0);
  const Scalar t3 = a(1) * a(1);
  const Scalar t4 = a(2) * a(2);
  return (t2 + t3 * 2.0 + t4 * 2.0 - a(0) * a(1) - a(0) * a(2) + sf * t2 -
          sf * t3 - sf * t4 + a(0) * a(1) * sf * 2.0 + a(0) * a(2) * sf * 2.0) /
         std::max(1e-63, t2 + t3 + t4);
}

inline Scalar get_rotated_pore_coeff_y(const Vec3x& a, const Scalar& sf) {
  const Scalar t2 = a(0) * a(0);
  const Scalar t3 = a(1) * a(1);
  const Scalar t4 = a(2) * a(2);
  return (t2 * 2.0 + t3 + t4 * 2.0 - a(0) * a(1) - a(1) * a(2) - sf * t2 +
          sf * t3 - sf * t4 + a(0) * a(1) * sf * 2.0 + a(1) * a(2) * sf * 2.0) /
         std::max(1e-63, t2 + t3 + t4);
}

inline Scalar get_rotated_pore_coeff_z(const Vec3x& a, const Scalar& sf) {
  const Scalar t2 = a(0) * a(0);
  const Scalar t3 = a(1) * a(1);
  const Scalar t4 = a(2) * a(2);
  return (t2 * 2.0 + t3 * 2.0 + t4 - a(0) * a(2) - a(1) * a(2) - sf * t2 -
          sf * t3 + sf * t4 + a(0) * a(2) * sf * 2.0 + a(1) * a(2) * sf * 2.0) /
         std::max(1e-63, t2 + t3 + t4);
}

inline Scalar get_rotated_drag_x(const Vec3x& a, const Vec3x& d) {
  const Scalar t2 = a(0) * a(0);
  const Scalar t3 = a(1) * a(1);
  const Scalar t4 = a(2) * a(2);
  const Scalar t5 = t2 + t3 + t4;
  const Scalar t8 = sqrt(t5);
  const Scalar t9 = t3 * t8;
  const Scalar t10 = a(2) * t2;
  const Scalar t6 = t9 + t10;
  const Scalar t7 = 1.0 / std::max(1e-63, t5);
  const Scalar t11 = t2 + t3;
  const Scalar t12 = a(2) - t8;
  const Scalar t13 = 1.0 / std::max(1e-63, t11 * t11);
  const Scalar t14 = 1.0 / std::max(1e-63, t11);
  return d(2) * t2 * t7 + d(0) * (t6 * t6) * t7 * t13 +
         a(0) * a(1) * d(2) * t7 + a(0) * a(2) * d(2) * t7 -
         a(0) * d(0) * t6 * t7 * t14 + d(1) * t2 * t3 * t7 * (t12 * t12) * t13 -
         a(0) * d(1) * t3 * t7 * t12 * t14 +
         a(0) * a(1) * d(0) * t6 * t7 * t12 * t13 +
         a(0) * a(1) * d(1) * t7 * t12 * t13 * (a(2) * t3 + t2 * t8);
}

inline Scalar get_rotated_drag_y(const Vec3x& a, const Vec3x& d) {
  const Scalar t2 = a(1) * a(1);
  const Scalar t3 = a(0) * a(0);
  const Scalar t4 = a(2) * a(2);
  const Scalar t5 = t2 + t3 + t4;
  const Scalar t8 = sqrt(t5);
  const Scalar t9 = t3 * t8;
  const Scalar t10 = a(2) * t2;
  const Scalar t6 = t9 + t10;
  const Scalar t7 = 1.0 / std::max(1e-63, t5);
  const Scalar t11 = t2 + t3;
  const Scalar t12 = a(2) - t8;
  const Scalar t13 = 1.0 / std::max(1e-63, t11 * t11);
  const Scalar t14 = 1.0 / std::max(1e-63, t11);
  return d(2) * t2 * t7 + d(1) * (t6 * t6) * t7 * t13 +
         a(0) * a(1) * d(2) * t7 + a(1) * a(2) * d(2) * t7 -
         a(1) * d(1) * t6 * t7 * t14 + d(0) * t2 * t3 * t7 * (t12 * t12) * t13 -
         a(1) * d(0) * t3 * t7 * t12 * t14 +
         a(0) * a(1) * d(1) * t6 * t7 * t12 * t13 +
         a(0) * a(1) * d(0) * t7 * t12 * t13 * (a(2) * t3 + t2 * t8);
}

inline Scalar get_rotated_drag_z(const Vec3x& a, const Vec3x& d) {
  const Scalar t2 = a(0) * a(0);
  const Scalar t3 = a(1) * a(1);
  const Scalar t4 = a(2) * a(2);
  const Scalar t5 = t2 + t3 + t4;
  const Scalar t6 = sqrt(t5);
  return (d(0) * (t2 * t2) + d(1) * (t3 * t3) + d(0) * t2 * t3 +
          d(1) * t2 * t3 + d(2) * t2 * t4 + d(2) * t3 * t4 -
          a(0) * d(0) * t3 * t6 + a(1) * d(0) * t2 * t6 +
          a(0) * d(1) * t3 * t6 - a(1) * d(1) * t2 * t6 -
          a(0) * a(2) * d(0) * t2 - a(1) * a(2) * d(0) * t2 -
          a(0) * a(2) * d(1) * t3 - a(1) * a(2) * d(1) * t3 +
          a(0) * a(2) * d(2) * t2 + a(0) * a(2) * d(2) * t3 +
          a(1) * a(2) * d(2) * t2 + a(1) * a(2) * d(2) * t3) /
         std::max(1e-63, t5 * (t2 + t3));
}

Scalar point_triangle_distance(const Vec3x& p, const Vec3x& a, const Vec3x& b,
                               const Vec3x& c, Scalar& t1, Scalar& t2,
                               Scalar& t3);

template <typename Callable, typename Array3T>
inline Scalar interpolate_2nd_order(const Vec3x& p, const Array3T& arr,
                                    Callable func) {
  Vec3x frac_p =
      Vec3x(p(0) - floor(p(0)), p(1) - floor(p(1)), p(2) - floor(p(2)));
  Vec3i base_p = Vec3i((int)floor(p(0)) + (frac_p(0) < 0.5 ? -1 : 0),
                       (int)floor(p(1)) + (frac_p(1) < 0.5 ? -1 : 0),
                       (int)floor(p(2)) + (frac_p(2) < 0.5 ? -1 : 0));

  Vec3x local_p =
      frac_p + Vec3x(frac_p(0) < 0.5 ? 1 : 0, frac_p(1) < 0.5 ? 1 : 0,
                     frac_p(2) < 0.5 ? 1 : 0);

  Scalar sum_u = 0.0;

  for (int k = 0; k < 3; ++k)
    for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 3; ++i) {
        Vec3i inode = base_p + Vec3i(i, j, k);
        if (!func(inode)) continue;

        Vec3x diff = local_p - Vec3x(i, j, k);

        sum_u += N_kernel<2>(diff) * arr(inode(0), inode(1), inode(2));
      }

  return sum_u;
}

// swap so that a<b
template <class T>
inline void sort(T& a, T& b) {
  if (a > b) std::swap(a, b);
}

// swap so that a<b<c
template <class T>
inline void sort(T& a, T& b, T& c) {
  if (a > b) std::swap(a, b);
  if (a > c) std::swap(a, c);
  if (b > c) std::swap(b, c);
}

// swap so that a<b<c<d
template <class T>
inline void sort(T& a, T& b, T& c, T& d) {
  if (a > b) std::swap(a, b);
  if (c > d) std::swap(c, d);
  if (a > c) std::swap(a, c);
  if (b > d) std::swap(b, d);
  if (b > c) std::swap(b, c);
}

inline void dhdr_yarn(Scalar mu, Scalar la, Scalar r22, Scalar r23, Scalar r33,
                      Scalar* dhdr22, Scalar* dhdr23, Scalar* dhdr33) {
  Eigen::Matrix2d r(2, 2);
  r(0, 0) = r22;
  r(0, 1) = r23;
  r(1, 0) = 0.0;
  r(1, 1) = r33;

  Eigen::JacobiSVD<Eigen::Matrix2d> svd(
      r, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Vector2d sigm = svd.singularValues();
  Eigen::Vector2d sigm_inv = Eigen::Vector2d(1.0 / sigm(0), 1.0 / sigm(1));
  Eigen::Vector2d lnsigm = Eigen::Vector2d(log(sigm(0)), log(sigm(1)));

  if (lnsigm(0) + lnsigm(1) > 0.0) {
    *dhdr22 = 0.0;
    *dhdr23 = 0.0;
    *dhdr33 = 0.0;
  } else {
    Eigen::Matrix2d tmp = Eigen::Matrix2d(sigm_inv.asDiagonal()) *
                          Eigen::Matrix2d(lnsigm.asDiagonal());
    Eigen::Matrix2d ret =
        svd.matrixU() *
        (2.0 * mu * tmp +
         la * lnsigm.sum() * Eigen::Matrix2d(sigm_inv.asDiagonal())) *
        svd.matrixV().transpose();

    *dhdr22 = ret(0, 0);
    *dhdr23 = ret(0, 1);
    *dhdr33 = ret(1, 1);
  }
}

inline void dgdr_cloth(Scalar mu, Scalar r13, Scalar r23, Scalar* dhdr13,
                       Scalar* dhdr23) {
  *dhdr13 = mu * r13;
  *dhdr23 = mu * r23;
}

inline void dhdr_cloth(Scalar mu, Scalar la, Scalar r33, Scalar* dhdr33) {
  if (r33 <= 1.0)
    *dhdr33 = -(2.0 * mu + la) * (1.0 - r33) * (1.0 - r33);
  else
    *dhdr33 = 0.0;
}

inline Scalar twist_component(const Eigen::Quaternion<Scalar>& rot,
                              const Vec3x& dir) {
  Vec3x ra = rot.vec();
  Vec3x p = dir * (ra.dot(dir));
  Eigen::Quaternion<Scalar> twist(rot.w(), p(0), p(1), p(2));
  twist.normalize();
  return Eigen::AngleAxis<Scalar>(twist).angle();
}

// build orthonormal basis from a unit vector
// Listing 2 in Frisvad 2012, Building an Orthonormal Basis from a 3D Unit
// Vector Without Normalization
inline void orthonormal_basis(const Vec3x& n, Vec3x& b1, Vec3x& b2) {
  if (n(2) + 1.0 < 1e-12) {
    b1 = Vec3x(0.0, -1.0, 0.0);
    b2 = Vec3x(-1.0, 0.0, 0.0);
    return;
  }

  const Scalar a = 1.0 / (1.0 + n(2));
  const Scalar b = -n(0) * n(1) * a;
  b1 = Vec3x(1.0 - n(0) * n(0) * a, b, -n(0));
  b2 = Vec3x(b, 1.0 - n(1) * n(1) * a, -n(1));
}

void print_histogram_analysis(const std::vector<VecXx>&, int,
                              const std::string&, bool);

template <typename T>
int bs_upper_bound(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a, T x) {
  int l = 0;
  int h = a.size();  // Not n - 1
  while (l < h) {
    int mid = (l + h) / 2;
    if (x >= a[mid]) {
      l = mid + 1;
    } else {
      h = mid;
    }
  }
  return l;
}

template <typename T>
T find_in_ordered(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, const T& x) {
  if (vec.size() == 0) return -1.0;

  int l = bs_upper_bound(vec, x) - 1;
  if (l < 0) return (T)l;
  if (l >= vec.size() - 1) {
    if (x > vec[vec.size() - 1])
      return vec.size();
    else
      return l;
  }

  const T d = vec[l + 1] - vec[l];
  if (d == (T)0.0) return (T)l + (T)0.5;

  return (T)l + (x - vec[l]) / d;
}

template <typename T>
T find_in_ordered(const std::vector<T>& vec, const T& x) {
  if (vec.size() == 0) return -1.0;

  int l = (int)(std::upper_bound(vec.begin(), vec.end(), x) - vec.begin()) - 1;
  if (l < 0) return (T)l;
  if (l >= vec.size() - 1) {
    if (x > vec[vec.size() - 1])
      return vec.size();
    else
      return l;
  }

  const T d = vec[l + 1] - vec[l];
  if (d == (T)0.0) return (T)l + (T)0.5;

  return (T)l + (x - vec[l]) / d;
}

template <typename T>
inline T cyl_h_from_area(T ra, T rb, T area) {
  const T pr = M_PI * (ra + rb);
  return (T)((sqrt(pr * pr + 4. * M_PI * area) - pr) / (2.0 * M_PI));
}

template <typename T>
inline T diffStrainElastic(T c0, T a, T h) {
  T t2 = c0 * (1.0 / 2.0);
  T t3 = c0 * c0;
  T t4 = t3 * (1.0 / 4.0);
  T t5 = t4 + 1.0;
  T t6 = sqrt(t5);
  T t7 = t2 + t6;
  return (exp(-a * h) * ((t7 * t7) * exp(a * h * 2.0) - 1.0)) / t7;
}

template <typename T, typename Callable>
inline T bisection_root_finding(T Ms, T ms, T eps, int max_iters, Callable g) {
  int count = 0;
  T s = (Ms + ms) * (T)0.5;

  while (1) {
    s = (Ms + ms) * (T)0.5;
    T v = g(s);

    if (count > max_iters || fabs(v) < eps * Ms) break;
    if (v > 0.0)
      Ms = s;
    else
      ms = s;

    count++;
  }

  return s;
}
}  // namespace strandsim

#endif
