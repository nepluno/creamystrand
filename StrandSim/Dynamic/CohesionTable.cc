/**
 * \copyright 2017 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "CohesionTable.hh"

#include <iostream>
#include <numeric>

#include "../Utils/MathUtilities.hh"

#define COLLISION_COHESION_TABLE

namespace strandsim {

inline Scalar spring_func(const Scalar& x, const Scalar& dmin, const Scalar& k0,
                          const Scalar& k1, const Scalar& x1, const Scalar& h) {
  const Scalar t2 = 1.0 / h;
  const Scalar t3 = x * x;
  const Scalar t4 = dmin + h;
  const Scalar t5 = h * h;
  const Scalar t6 = k1 * t5;
  const Scalar t7 = dmin * h * k1;
  return -t2 - t2 * t3 * 1.0 / (t4 * t4) * (t6 + t7 - h * x1 * 3.0 - 3.0) +
         t2 * t3 * 1.0 / (t4 * t4 * t4) * x * (t6 + t7 - h * x1 * 2.0 - 2.0);
}

inline Scalar grad_spring_func(const Scalar& x, const Scalar& dmin,
                               const Scalar& k0, const Scalar& k1,
                               const Scalar& x1, const Scalar& h) {
  const Scalar t2 = 1.0 / h;
  const Scalar t3 = dmin + h;
  const Scalar t4 = h * h;
  const Scalar t5 = k1 * t4;
  const Scalar t6 = dmin * h * k1;
  return t2 * 1.0 / (t3 * t3 * t3) * (x * x) * (t5 + t6 - h * x1 * 2.0 - 2.0) *
             3.0 -
         t2 * 1.0 / (t3 * t3) * x * (t5 + t6 - h * x1 * 3.0 - 3.0) * 2.0;
}

CohesionTable::CohesionTable()
    : m_sigma(72.0),
      m_theta(M_PI / 4),
      m_radii(0.004),
      m_max_alpha(M_PI - m_theta),
      m_max_d0(0.2),
      m_min_d0(2.0 * m_radii),
      m_min_d0_planar(m_radii),
      m_discretization(256) {}

void CohesionTable::print_table(std::ostream& oss, const MatXx& mat,
                                const Scalar& dmin) const {
  Scalar d_inc = m_max_d0 / (Scalar)m_discretization;

  for (int j = 0; j < m_discretization; ++j) {
    for (int i = 0; i < m_discretization; ++i) {
      Scalar d0 = (Scalar)i * d_inc + dmin;
      Scalar a_inc = (m_max_As(i) - m_min_As(i)) / m_discretization;

      Scalar A = (Scalar)j * a_inc + m_min_As(i);

      Scalar value = mat(j, i);

      oss << A << " " << d0 << " " << value << std::endl;
    }
  }
}

void CohesionTable::print_dEdd_table(std::ostream& oss) const {
  print_table(oss, m_dEdd_table, 0.0);
}

Scalar CohesionTable::interpolate_table_planar(const Scalar& A,
                                               const Scalar& d0,
                                               const MatXx& mat,
                                               const Scalar& dmin) const {
  Scalar d_inc = m_max_d0 / (Scalar)m_discretization;
  Scalar p = std::max(0., d0 - dmin) / d_inc;
  Scalar fp = std::min(p - floor(p), 1.0);
  int ip0 = std::max(0, std::min(m_discretization - 2, (int)p));
  int ip1 = ip0 + 1;

  Scalar a_inc0 =
      (m_max_A_planars(ip0) - m_min_A_planars(ip0)) / m_discretization;
  Scalar a_inc1 =
      (m_max_A_planars(ip1) - m_min_A_planars(ip1)) / m_discretization;

  if (a_inc0 == 0.0 || a_inc1 == 0.0) return 0.0;

  Scalar q0 =
      (std::max(m_min_A_planars(ip0), std::min(m_max_A_planars(ip0), A)) -
       m_min_A_planars(ip0)) /
      a_inc0;
  Scalar q1 =
      (std::max(m_min_A_planars(ip1), std::min(m_max_A_planars(ip1), A)) -
       m_min_A_planars(ip1)) /
      a_inc1;

  Scalar fq0 = std::min(q0 - floor(q0), 1.0);
  Scalar fq1 = std::min(q1 - floor(q1), 1.0);

  int iq00 = std::max(0, std::min(m_discretization - 2, (int)q0));
  int iq10 = std::max(0, std::min(m_discretization - 2, (int)q1));

  int iq01 = iq00 + 1;
  int iq11 = iq10 + 1;

  Scalar dEdd0;
  Scalar dEdd1;

  Scalar v00 = mat(iq00, ip0);
  Scalar v01 = mat(iq01, ip0);
  Scalar v10 = mat(iq10, ip1);
  Scalar v11 = mat(iq11, ip1);

  dEdd0 = lerp(v00, v01, fq0);
  dEdd1 = lerp(v10, v11, fq1);

  return lerp(dEdd0, dEdd1, fp);
}

Scalar CohesionTable::interpolate_table(const Scalar& A, const Scalar& d0,
                                        const MatXx& mat,
                                        const Scalar& dmin) const {
  Scalar d_inc = m_max_d0 / (Scalar)m_discretization;
  Scalar p = std::max(0., d0 - dmin) / d_inc;
  Scalar fp = std::min(p - floor(p), 1.0);
  int ip0 = std::max(0, std::min(m_discretization - 2, (int)p));
  int ip1 = ip0 + 1;

  Scalar a_inc0 = (m_max_As(ip0) - m_min_As(ip0)) / m_discretization;
  Scalar a_inc1 = (m_max_As(ip1) - m_min_As(ip1)) / m_discretization;

  if (a_inc0 == 0.0 || a_inc1 == 0.0) return 0.0;

  Scalar q0 =
      (std::max(m_min_As(ip0), std::min(m_max_As(ip0), A)) - m_min_As(ip0)) /
      a_inc0;
  Scalar q1 =
      (std::max(m_min_As(ip1), std::min(m_max_As(ip1), A)) - m_min_As(ip1)) /
      a_inc1;

  Scalar fq0 = std::min(q0 - floor(q0), 1.0);
  Scalar fq1 = std::min(q1 - floor(q1), 1.0);

  int iq00 = std::max(0, std::min(m_discretization - 2, (int)q0));
  int iq10 = std::max(0, std::min(m_discretization - 2, (int)q1));

  int iq01 = iq00 + 1;
  int iq11 = iq10 + 1;

  Scalar dEdd0;
  Scalar dEdd1;

  Scalar v00 = mat(iq00, ip0);
  Scalar v01 = mat(iq01, ip0);
  Scalar v10 = mat(iq10, ip1);
  Scalar v11 = mat(iq11, ip1);

  dEdd0 = lerp(v00, v01, fq0);
  dEdd1 = lerp(v10, v11, fq1);
  //
  //        if(std::isnan(dEdd0)) {
  //            std::cout << "dEdd0 is NAN!" << std::endl;
  //            std::cout << "A: " << A << std::endl;
  //            std::cout << "d0: " << d0 << std::endl;
  //            std::cout << "mat: " << mat << std::endl;
  //            std::cout << "dmin: " << dmin << std::endl;
  //            std::cout << "fp: " << fp << std::endl;
  //
  //            exit(0);
  //        }
  //
  //        if(std::isnan(dEdd1)) {
  //            std::cout << "dEdd1 is NAN!" << std::endl;
  //            std::cout << "A: " << A << std::endl;
  //            std::cout << "d0: " << d0 << std::endl;
  //            std::cout << "mat: " << mat << std::endl;
  //            std::cout << "dmin: " << dmin << std::endl;
  //            std::cout << "fp: " << fp << std::endl;
  //
  //
  //            exit(0);
  //        }
  //
  return lerp(dEdd0, dEdd1, fp);
}

Scalar CohesionTable::interpolate_table_grad_planar(const Scalar& A,
                                                    const Scalar& d0,
                                                    const MatXx& mat,
                                                    const Scalar& dmin) const {
  Scalar d_inc = m_max_d0 / (Scalar)m_discretization;
  Scalar p = std::max(0., d0 - dmin) / d_inc;
  int ip0 = std::max(0, std::min(m_discretization - 2, (int)p));
  int ip1 = ip0 + 1;

  Scalar a_inc0 =
      (m_max_A_planars(ip0) - m_min_A_planars(ip0)) / m_discretization;
  Scalar a_inc1 =
      (m_max_A_planars(ip1) - m_min_A_planars(ip1)) / m_discretization;

  if (a_inc0 == 0.0 || a_inc1 == 0.0) return 0.0;

  Scalar q0 =
      (std::max(m_min_A_planars(ip0), std::min(m_max_A_planars(ip0), A)) -
       m_min_A_planars(ip0)) /
      a_inc0;
  Scalar q1 =
      (std::max(m_min_A_planars(ip1), std::min(m_max_A_planars(ip1), A)) -
       m_min_A_planars(ip1)) /
      a_inc1;

  Scalar fq0 = std::min(q0 - floor(q0), 1.0);
  Scalar fq1 = std::min(q1 - floor(q1), 1.0);

  int iq00 = std::max(0, std::min(m_discretization - 2, (int)q0));
  int iq10 = std::max(0, std::min(m_discretization - 2, (int)q1));

  int iq01 = iq00 + 1;
  int iq11 = iq10 + 1;

  Scalar dEdd0;
  Scalar dEdd1;

  Scalar v00 = mat(iq00, ip0);
  Scalar v01 = mat(iq01, ip0);
  Scalar v10 = mat(iq10, ip1);
  Scalar v11 = mat(iq11, ip1);

  dEdd0 = lerp(v00, v01, fq0);
  dEdd1 = lerp(v10, v11, fq1);

  return (dEdd1 - dEdd0) / d_inc;
}

Scalar CohesionTable::interpolate_table_grad(const Scalar& A, const Scalar& d0,
                                             const MatXx& mat,
                                             const Scalar& dmin) const {
  Scalar d_inc = m_max_d0 / (Scalar)m_discretization;
  Scalar p = std::max(0., d0 - dmin) / d_inc;
  int ip0 = std::max(0, std::min(m_discretization - 2, (int)p));
  int ip1 = ip0 + 1;

  Scalar a_inc0 = (m_max_As(ip0) - m_min_As(ip0)) / m_discretization;
  Scalar a_inc1 = (m_max_As(ip1) - m_min_As(ip1)) / m_discretization;

  if (a_inc0 == 0.0 || a_inc1 == 0.0) return 0.0;

  Scalar q0 =
      (std::max(m_min_As(ip0), std::min(m_max_As(ip0), A)) - m_min_As(ip0)) /
      a_inc0;
  Scalar q1 =
      (std::max(m_min_As(ip1), std::min(m_max_As(ip1), A)) - m_min_As(ip1)) /
      a_inc1;

  Scalar fq0 = std::min(q0 - floor(q0), 1.0);
  Scalar fq1 = std::min(q1 - floor(q1), 1.0);

  int iq00 = std::max(0, std::min(m_discretization - 2, (int)q0));
  int iq10 = std::max(0, std::min(m_discretization - 2, (int)q1));

  int iq01 = iq00 + 1;
  int iq11 = iq10 + 1;

  Scalar dEdd0;
  Scalar dEdd1;

  Scalar v00 = mat(iq00, ip0);
  Scalar v01 = mat(iq01, ip0);
  Scalar v10 = mat(iq10, ip1);
  Scalar v11 = mat(iq11, ip1);

  dEdd0 = lerp(v00, v01, fq0);
  dEdd1 = lerp(v10, v11, fq1);

  return (dEdd1 - dEdd0) / d_inc;
}

Scalar CohesionTable::interpolate_dEdd(const Scalar& A,
                                       const Scalar& d0) const {
  return interpolate_table(A, d0, m_dEdd_table, 0.0);
}

Scalar CohesionTable::interpolate_d2Edd2(const Scalar& A,
                                         const Scalar& d0) const {
  return interpolate_table_grad(A, d0, m_dEdd_table, 0.0);
}

Scalar CohesionTable::interpolate_alpha(const Scalar& A,
                                        const Scalar& d0) const {
  return interpolate_table(A, d0, m_alpha_table, m_min_d0);
}

Scalar CohesionTable::interpolate_dEdd_planar(const Scalar& A,
                                              const Scalar& d0) const {
  return interpolate_table_planar(A, d0, m_dEdd_planar_table, 0.0);
}

Scalar CohesionTable::interpolate_d2Edd2_planar(const Scalar& A,
                                                const Scalar& d0) const {
  return interpolate_table_grad_planar(A, d0, m_dEdd_planar_table, 0.0);
}

Scalar CohesionTable::interpolate_alpha_planar(const Scalar& A,
                                               const Scalar& d0) const {
  return interpolate_table_planar(A, d0, m_alpha_planar_table, m_min_d0_planar);
}

void CohesionTable::setParameter(const Scalar& sigma, const Scalar& theta,
                                 const Scalar& radii, const Scalar& max_d0,
                                 const int disc) {
  m_sigma = sigma;
  m_theta = theta;
  m_radii = radii;
  m_max_d0 = max_d0;
  m_discretization = disc;
  m_min_d0 = radii * 2.0;
  m_max_alpha = M_PI - m_theta;
}

Scalar CohesionTable::getRadii() const { return m_radii; }

Scalar CohesionTable::computeH(const Scalar& R, const Scalar& alpha) const {
  return m_radii * sin(alpha) - R * (1.0 - sin(m_theta + alpha));
}

Scalar CohesionTable::computeR(const Scalar& alpha, const Scalar& d0) const {
  return (d0 - 2.0 * m_radii * cos(alpha)) / (2.0 * cos(m_theta + alpha));
}

Scalar CohesionTable::computeA(const Scalar& R, const Scalar& alpha) const {
  return 2.0 * R * R *
             (alpha + m_theta - M_PI / 2 + 0.5 * sin(2.0 * (alpha + m_theta))) +
         2.0 * m_radii * R * (sin(2.0 * alpha + m_theta) - sin(m_theta)) -
         m_radii * m_radii * (2.0 * alpha - sin(2.0 * alpha));
}

Scalar CohesionTable::computeApproxA(const Scalar& alpha,
                                     const Scalar& d0) const {
  const Scalar gamma = alpha + m_theta;
  const Scalar t2 = m_radii * m_radii;
  const Scalar t3 = d0 * d0;
  const Scalar t4 = cos(m_theta);
  const Scalar t5 = t4 * t4;
  const Scalar t6 = sin(m_theta);
  return -t2 * sin(m_theta * 2.0) + t2 * M_PI * (1.0 / 3.0) -
         t3 * M_PI * (1.0 / 6.0) - gamma * t2 * (8.0 / 3.0) +
         gamma * t3 * (1.0 / 3.0) + t2 * m_theta * 2.0 -
         t2 * t5 * M_PI * (4.0 / 3.0) + d0 * m_radii * t4 * 2.0 +
         gamma * t2 * t5 * (8.0 / 3.0) +
         d0 * gamma * m_radii * t6 * (2.0 / 3.0) -
         d0 * m_radii * t6 * M_PI * (1.0 / 3.0);
}

Scalar CohesionTable::computeApproxdEdd(const Scalar& alpha,
                                        const Scalar& d0) const {
  const Scalar gamma = alpha + m_theta;

  const Scalar t2 = sin(m_theta);
  const Scalar t3 = m_radii * m_radii;
  const Scalar t4 = cos(m_theta);
  const Scalar t5 = d0 * d0;
  const Scalar t6 = d0 * m_radii * t2 * 2.0;
  const Scalar t7 = m_theta * 2.0;
  const Scalar t8 = sin(t7);
  return (m_sigma *
          (t3 * -8.0 + t5 + t6 + t3 * (t4 * t4) * 8.0 + t3 * t8 * M_PI * 2.0 -
           gamma * t3 * t8 * 4.0 - d0 * gamma * m_radii * t4 * 2.0 +
           d0 * m_radii * t4 * M_PI) *
          2.0) /
         (t5 + t6 - (t2 * t2) * t3 * 8.0);
}

Scalar CohesionTable::computedEddPlanar(const Scalar& R,
                                        const Scalar& alpha) const {
  if (R == 0.0) {
    return 0.0;
  } else {
    const Scalar t2 = m_theta * 2.0;
    const Scalar t3 = sin(alpha);
    const Scalar t4 = alpha + m_theta;
    const Scalar t5 = sin(t4);
    const Scalar t6 = alpha + t2;
    const Scalar t7 = m_radii * m_radii;
    const Scalar t8 = sin(t6);
    const Scalar t9 = R * R;
    const Scalar t10 = alpha * 2.0;
    const Scalar t11 = sin(t2);
    const Scalar t12 = t2 + t10;
    const Scalar t13 = sin(t12);
    const Scalar t14 = m_theta * 3.0;
    const Scalar t15 = alpha + t14;
    const Scalar t16 = cos(t10);
    const Scalar t17 = cos(t12);
    const Scalar t18 = cos(m_theta);
    const Scalar t19 = t10 + m_theta;
    const Scalar t20 = cos(t19);
    return (m_sigma *
            (t7 * M_PI * -2.0 - t9 * M_PI * 2.0 + alpha * t7 * 2.0 +
             alpha * t9 * 2.0 + t3 * t7 * 4.0 + t3 * t9 * 4.0 + t7 * t8 * 2.0 +
             t8 * t9 * 4.0 - t7 * t11 * 2.0 + t7 * t13 * 2.0 - t9 * t11 +
             t9 * t13 * 3.0 + t7 * m_theta * 4.0 + t9 * m_theta * 4.0 +
             t7 * sin(t10) * 2.0 + t7 * sin(alpha - t2) * 2.0 +
             t9 * sin(t10 + m_theta * 4.0) - R * m_radii * sin(t14) +
             R * m_radii * sin(t15) * 2.0 + R * m_radii * sin(t19) * 5.0 -
             R * m_radii * sin(m_theta) * 3.0 +
             R * m_radii * sin(alpha - m_theta) * 6.0 + t7 * t16 * M_PI * 2.0 +
             t9 * t17 * M_PI * 2.0 + R * m_radii * t5 * 8.0 -
             alpha * t7 * t16 * 2.0 - alpha * t9 * t17 * 2.0 -
             t7 * t16 * m_theta * 4.0 - t9 * t17 * m_theta * 4.0 +
             R * m_radii * sin(t10 + t14) * 3.0 +
             R * m_radii * t18 * m_theta * 8.0 -
             R * m_radii * t20 * m_theta * 8.0 -
             R * m_radii * t18 * M_PI * 4.0 + R * m_radii * t20 * M_PI * 4.0 +
             R * alpha * m_radii * t18 * 4.0 -
             R * alpha * m_radii * t20 * 4.0)) /
           (R * (m_radii * 2.0 + R * t18 * 4.0 + R * cos(t4) * 3.0 +
                 R * cos(t15) + m_radii * cos(alpha) * 2.0 +
                 m_radii * cos(t2) * 2.0 + m_radii * cos(t6) * 2.0 -
                 R * t5 * M_PI * 2.0 + R * alpha * t5 * 2.0 -
                 m_radii * t3 * M_PI * 2.0 + R * t5 * m_theta * 4.0 +
                 alpha * m_radii * t3 * 2.0 + m_radii * t3 * m_theta * 4.0));
  }
}

Scalar CohesionTable::computedEdd(const Scalar& R, const Scalar& alpha) const {
  if (R == 0.0) {
    return 0.0;
  } else {
    const Scalar t2 = sin(alpha);
    const Scalar t3 = alpha + m_theta;
    const Scalar t4 = sin(t3);
    const Scalar t5 = R * R;
    const Scalar t6 = m_radii * m_radii;
    const Scalar t7 = m_theta * 2.0;
    const Scalar t8 = alpha * 2.0;
    const Scalar t9 = t7 + t8;
    const Scalar t10 = sin(t9);
    const Scalar t11 = cos(t8);
    const Scalar t12 = cos(t9);
    const Scalar t13 = cos(m_theta);
    const Scalar t14 = t8 + m_theta;
    const Scalar t15 = cos(t14);
    return (m_sigma *
            (-t5 * M_PI - t6 * M_PI + alpha * t5 * 2.0 + alpha * t6 * 2.0 +
             t5 * t10 * 2.0 + t6 * t10 + t5 * m_theta * 2.0 +
             t6 * m_theta * 2.0 - t6 * sin(t7) + t6 * sin(t8) +
             R * m_radii * sin(t14) * 3.0 - R * m_radii * sin(m_theta) * 2.0 +
             R * m_radii * sin(t8 + m_theta * 3.0) + t5 * t12 * M_PI +
             t6 * t11 * M_PI - alpha * t5 * t12 * 2.0 - alpha * t6 * t11 * 2.0 -
             t5 * t12 * m_theta * 2.0 - t6 * t11 * m_theta * 2.0 +
             R * m_radii * t13 * m_theta * 4.0 -
             R * m_radii * t15 * m_theta * 4.0 -
             R * m_radii * t13 * M_PI * 2.0 + R * m_radii * t15 * M_PI * 2.0 +
             R * alpha * m_radii * t13 * 4.0 -
             R * alpha * m_radii * t15 * 4.0)) /
           (R * (m_radii * cos(alpha + t7) + R * cos(t3) * 2.0 +
                 m_radii * cos(alpha) - R * t4 * M_PI + R * alpha * t4 * 2.0 -
                 m_radii * t2 * M_PI + R * t4 * m_theta * 2.0 +
                 alpha * m_radii * t2 * 2.0 + m_radii * t2 * m_theta * 2.0));
  }
}

Scalar CohesionTable::computeRPlanar(const Scalar& alpha,
                                     const Scalar& d0) const {
  return (d0 - m_radii * cos(alpha)) / (cos(m_theta + alpha) + cos(m_theta));
}

Scalar CohesionTable::computeAPlanar(const Scalar& R,
                                     const Scalar& alpha) const {
  return 2.0 * (0.5 * m_radii * m_radii * sin(alpha) * cos(alpha) +
                m_radii * sin(alpha) * R * cos(m_theta + alpha) +
                0.5 * R * R * sin(m_theta + alpha) * cos(m_theta + alpha)) +
         2.0 *
             (R * sin(m_theta + alpha) - R * sin(m_theta) +
              m_radii * sin(alpha)) *
             R * cos(m_theta) +
         R * R * sin(m_theta) * cos(m_theta) -
         (alpha * m_radii * m_radii + R * R * (M_PI - 2.0 * m_theta - alpha));
}

Scalar CohesionTable::computeHPlanar(const Scalar& R,
                                     const Scalar& alpha) const {
  return computeH(R, alpha);
}

Scalar CohesionTable::computeApproxAPlanar(const Scalar& alpha,
                                           const Scalar& d0) const {
  if (m_theta == 0.0)
    return 0.0;
  else {
    const Scalar gamma = alpha + m_theta * 2.0;
    const Scalar t2 = m_theta * 3.0;
    const Scalar t3 = cos(t2);
    const Scalar t4 = m_radii * m_radii;
    const Scalar t5 = d0 * d0;
    const Scalar t6 = cos(m_theta);
    const Scalar t7 = M_PI * M_PI;
    const Scalar t8 = m_theta * 5.0;
    const Scalar t9 = cos(t8);
    const Scalar t10 = gamma * gamma;
    const Scalar t11 = sin(m_theta);
    const Scalar t12 = sin(t2);
    const Scalar t13 = sin(t8);
    return 1.0 / (t11 * t11 * t11) *
           (t3 * t4 * 6.0 - t3 * t5 * 1.2E1 + t5 * t6 * 1.2E1 - t4 * t9 * 6.0 -
            t4 * t11 * M_PI * 4.0 - t4 * t12 * M_PI * 1.6E1 -
            t5 * t11 * M_PI * 3.2E1 + t4 * t13 * M_PI * 4.0 -
            d0 * m_radii * t3 * 2.4E1 + d0 * m_radii * t6 * 2.4E1 -
            gamma * t4 * t11 * 3.2E1 + gamma * t4 * t12 * 2.8E1 +
            gamma * t5 * t11 * 3.2E1 - gamma * t4 * t13 * 4.0 -
            t3 * t4 * t7 * 3.0 - t3 * t4 * t10 * 3.0 + t4 * t6 * t7 * 2.2E1 +
            t5 * t6 * t7 * 2.0E1 + t4 * t6 * t10 * 2.2E1 + t4 * t7 * t9 +
            t5 * t6 * t10 * 2.0E1 + t4 * t9 * t10 + t4 * t11 * m_theta * 7.2E1 -
            t4 * t12 * m_theta * 2.4E1 + d0 * gamma * m_radii * t11 * 4.0E1 +
            d0 * gamma * m_radii * t12 * 8.0 + d0 * m_radii * t6 * t7 * 4.0E1 +
            d0 * m_radii * t6 * t10 * 4.0E1 -
            d0 * m_radii * t11 * M_PI * 4.0E1 -
            d0 * m_radii * t12 * M_PI * 8.0 + gamma * t3 * t4 * M_PI * 6.0 -
            gamma * t4 * t6 * M_PI * 4.4E1 - gamma * t5 * t6 * M_PI * 4.0E1 -
            gamma * t4 * t9 * M_PI * 2.0 -
            d0 * gamma * m_radii * t6 * M_PI * 8.0E1) *
           (1.0 / 4.8E1);
  }
}

Scalar CohesionTable::computeApproxdEddPlanar(const Scalar& alpha,
                                              const Scalar& d0) const {
  const Scalar gamma = alpha + m_theta * 2.0;

  if (m_theta == 0.0) {
    const Scalar t2 = d0 * d0;
    const Scalar t3 = m_radii * m_radii;
    return m_sigma * 1.0 / pow(d0 + m_radii, 2.0) *
           (t2 * M_PI * 4.0 + t3 * M_PI * 4.0 - gamma * t2 * 4.0 -
            gamma * t3 * 4.0 + d0 * m_radii * M_PI * 8.0 -
            d0 * gamma * m_radii * 8.0) *
           (1.0 / 2.0);
  } else {
    const Scalar t2 = m_theta * 2.0;
    const Scalar t3 = cos(t2);
    const Scalar t4 = d0 * d0;
    const Scalar t5 = m_radii * m_radii;
    const Scalar t6 = t3 * t3;
    const Scalar t7 = M_PI * M_PI;
    const Scalar t8 = gamma * gamma;
    const Scalar t9 = sin(t2);
    const Scalar t10 = m_theta * 4.0;
    const Scalar t11 = sin(t10);
    return (m_sigma * 1.0 / pow(d0 + m_radii * t3, 2.0) *
            (t4 * -2.0 + t3 * t4 * 2.0 + t4 * t7 - t5 * t6 * 2.0 + t4 * t8 -
             t5 * t7 * 2.0 - t5 * t8 * 2.0 - gamma * t4 * M_PI * 2.0 +
             gamma * t5 * M_PI * 4.0 - t4 * t9 * M_PI * 2.0 - t5 * t11 * M_PI -
             d0 * m_radii * t3 * 4.0 + d0 * m_radii * t6 * 4.0 -
             d0 * m_radii * t7 - d0 * m_radii * t8 + gamma * t4 * t9 * 2.0 +
             gamma * t5 * t11 - t3 * t4 * t7 + t3 * t5 * t6 * 2.0 -
             t3 * t4 * t8 + t3 * t5 * t7 + t3 * t5 * t8 + t5 * t6 * t7 +
             t5 * t6 * t8 + d0 * gamma * m_radii * t9 * 2.0 +
             d0 * gamma * m_radii * t11 + d0 * m_radii * t6 * t7 +
             d0 * m_radii * t6 * t8 + d0 * gamma * m_radii * M_PI * 2.0 -
             d0 * m_radii * t9 * M_PI * 2.0 - d0 * m_radii * t11 * M_PI +
             gamma * t3 * t4 * M_PI * 2.0 - gamma * t3 * t5 * M_PI * 2.0 -
             gamma * t5 * t6 * M_PI * 2.0 -
             d0 * gamma * m_radii * t6 * M_PI * 2.0) *
            (-1.0 / 2.0)) /
           sin(m_theta);
  }
}

Scalar CohesionTable::computedEddAreaDist(const Scalar& A_target,
                                          const Scalar& d0) const {
  if (d0 < getDMin()) return computedEddAreaDist(A_target, getDMin());

  //        if(d0 < 2.0 * sqrt(A_target / M_PI + 2.0 * m_radii * m_radii) -
  //        m_radii * 2.0)
  //            return 0.0;

  if (d0 > (1.0 + 0.5 * m_theta) * sqrt(A_target) + m_radii * 2.0) return 0.0;

  Scalar alpha = interpolate_alpha(A_target, d0);
  Scalar gamma = alpha + m_theta;

  Scalar dEdd;

  if (gamma < M_PI / 2. + m_ang_epsilon && gamma > M_PI / 2. - m_ang_epsilon) {
    dEdd = computeApproxdEdd(alpha, d0);
  } else {
    Scalar R_target = computeR(alpha, d0);
    dEdd = computedEdd(R_target, alpha);
  }

  return std::max(0.0, dEdd);
}

Scalar CohesionTable::computedEddAreaDistPlanar(const Scalar& A_target,
                                                const Scalar& d0) const {
  if (d0 < getDMinPlanar())
    return computedEddAreaDistPlanar(A_target, getDMinPlanar());

  if (m_theta > 0.0 && d0 < sqrt((1.0 - cos(m_theta)) * (1.0 - cos(m_theta)) *
                                 (A_target + M_PI * m_radii * m_radii) /
                                 (m_theta - 0.5 * sin(2.0 * m_theta))) -
                                m_radii)
    return 0.0;

  Scalar alpha = interpolate_alpha_planar(A_target, d0);

  Scalar gamma = alpha + m_theta * 2.0;

  Scalar dEdd;

  if (gamma < M_PI + m_ang_epsilon && gamma > M_PI - m_ang_epsilon) {
    dEdd = computeApproxdEddPlanar(alpha, d0);
  } else {
    Scalar R_target = computeRPlanar(alpha, d0);
    dEdd = computedEddPlanar(R_target, alpha);
  }

  return std::max(0.0, dEdd);
}

void CohesionTable::construct_alpha_table() {
  Scalar alpha_inc = m_max_alpha / (Scalar)m_discretization;
  Scalar d_inc = m_max_d0 / (Scalar)m_discretization;

  m_alpha_table.resize(m_discretization, m_discretization);
  m_A_table.resize(m_discretization, m_discretization);
  m_dEdd_table.resize(m_discretization, m_discretization);
  m_max_As.resize(m_discretization);
  m_min_As.resize(m_discretization);

  const Scalar dmin = getDMin();

  for (int i = 0; i < m_discretization; ++i) {
    Scalar d0 = d_inc * (Scalar)i + dmin;

    Scalar maxA = -1e+20;
    Scalar minA = 1e+20;
    for (int j = 0; j < m_discretization; ++j) {
      Scalar alpha = alpha_inc * (Scalar)j;

      Scalar gamma = alpha + m_theta;

      Scalar A;

      if (gamma < M_PI / 2. + m_ang_epsilon &&
          gamma > M_PI / 2. - m_ang_epsilon) {
        A = computeApproxA(alpha, d0);
      } else {
        Scalar R = computeR(alpha, d0);
        A = computeA(R, alpha);
      }

      m_A_table(j, i) = A;
      maxA = std::max(maxA, std::max(0.0, A));
      minA = std::min(minA, std::max(0.0, A));
    }

    m_max_As(i) = maxA;
    m_min_As(i) = minA;
  }

  // std::cout << m_A_table << std::endl;

  for (int i = 0; i < m_discretization; ++i) {
    const VecXx& v = m_A_table.col(i);

    Scalar A_inc = (m_max_As(i) - m_min_As(i)) / m_discretization;

    for (int j = 0; j < m_discretization; ++j) {
      Scalar A_target = A_inc * (Scalar)j + m_min_As(i);

      Scalar alpha_target =
          clamp(find_in_ordered(v, A_target) * alpha_inc, 0.0, m_max_alpha);

      m_alpha_table(j, i) = alpha_target;
    }
  }

  //        MatXx buffer(m_discretization * m_discretization, 3);

  for (int j = 0; j < m_discretization; ++j) {
    for (int i = 0; i < m_discretization; ++i) {
      Scalar A_inc = (m_max_As(i) - m_min_As(i)) / m_discretization;
      Scalar A_target = std::max(m_min_As(i), A_inc * (Scalar)j + m_min_As(i));

      Scalar d0 = d_inc * (Scalar)i;

      m_dEdd_table(j, i) = computedEddAreaDist(A_target, d0);
      //
      //                buffer(j * m_discretization + i, 0) = A_target;
      //                buffer(j * m_discretization + i, 1) = d0;
      //                buffer(j * m_discretization + i, 2) = m_dEdd_table(j,
      //                i);
    }
  }

  //        std::cout << buffer << std::endl;

  //        std::cout << m_dEdd_table << std::endl;
  // std::cout << m_d2Edd2_table << std::endl;
}

void CohesionTable::construct_planar_alpha_table() {
  Scalar alpha_inc = m_max_alpha / (Scalar)m_discretization;
  Scalar d_inc = m_max_d0 / (Scalar)m_discretization;

  m_alpha_planar_table.resize(m_discretization, m_discretization);
  m_A_planar_table.resize(m_discretization, m_discretization);
  m_dEdd_planar_table.resize(m_discretization, m_discretization);
  m_max_A_planars.resize(m_discretization);
  m_min_A_planars.resize(m_discretization);

  const Scalar dmin = getDMinPlanar();

  for (int i = 0; i < m_discretization; ++i) {
    Scalar d0 = d_inc * (Scalar)i + dmin;

    Scalar maxA = -1e+20;
    Scalar minA = 1e+20;
    for (int j = 0; j < m_discretization; ++j) {
      Scalar alpha = alpha_inc * (Scalar)j;

      Scalar gamma = alpha + m_theta * 2.0;

      Scalar A;

      if (gamma < M_PI + m_ang_epsilon && gamma > M_PI - m_ang_epsilon) {
        A = computeApproxAPlanar(alpha, d0);
      } else {
        Scalar R = computeRPlanar(alpha, d0);
        A = computeAPlanar(R, alpha);
      }

      m_A_planar_table(j, i) = A;
      maxA = std::max(maxA, std::max(0.0, A));
      minA = std::min(minA, std::max(0.0, A));

      // std::cout << maxA << ", " << minA << ", " << alpha << ", " << j <<
      // std::endl;
    }

    m_max_A_planars(i) = maxA;
    m_min_A_planars(i) = minA;
  }

  // std::cout << m_A_planar_table << std::endl;

  for (int i = 0; i < m_discretization; ++i) {
    const VecXx& v = m_A_planar_table.col(i);
    const size_t num_results = 1;

    Scalar A_inc = (m_max_A_planars(i) - m_min_A_planars(i)) / m_discretization;

    for (int j = 0; j < m_discretization; ++j) {
      Scalar A_target = A_inc * (Scalar)j + m_min_A_planars(i);

      Scalar alpha_target =
          clamp(find_in_ordered(v, A_target) * alpha_inc, 0.0, m_max_alpha);

      m_alpha_planar_table(j, i) = alpha_target;
    }
  }

  for (int j = 0; j < m_discretization; ++j) {
    for (int i = 0; i < m_discretization; ++i) {
      Scalar A_inc =
          (m_max_A_planars(i) - m_min_A_planars(i)) / m_discretization;
      Scalar A_target =
          std::max(m_min_A_planars(i), A_inc * (Scalar)j + m_min_A_planars(i));

      Scalar d0 = d_inc * (Scalar)i;
      m_dEdd_planar_table(j, i) = computedEddAreaDistPlanar(A_target, d0);
    }
  }

  // std::cout << m_dEdd_table << std::endl;
  // std::cout << m_d2Edd2_table << std::endl;
}

Scalar CohesionTable::getDMinPlanar() const { return m_min_d0_planar; }

Scalar CohesionTable::getDMin() const { return m_min_d0; }

};  // namespace strandsim
