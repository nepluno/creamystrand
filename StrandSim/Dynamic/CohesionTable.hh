/**
 * \copyright 2017 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef __COHESION_TABLE_GEN_H__
#define __COHESION_TABLE_GEN_H__

#include <Eigen/Core>

#include "../Core/Definitions.hh"

namespace strandsim {
class CohesionTable {
  // const parameters
  Scalar m_sigma;
  Scalar m_theta;
  Scalar m_radii;
  Scalar m_max_alpha;
  Scalar m_max_d0;

  Scalar m_min_d0;
  Scalar m_min_d0_planar;

  const Scalar m_ang_epsilon = 0.008;

  int m_discretization;

  MatXx m_A_table;
  MatXx m_alpha_table;
  MatXx m_dEdd_table;

  MatXx m_A_planar_table;
  MatXx m_alpha_planar_table;
  MatXx m_dEdd_planar_table;

  VecXx m_max_As;
  VecXx m_min_As;

  VecXx m_max_A_planars;
  VecXx m_min_A_planars;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Scalar computeR(const Scalar& alpha, const Scalar& d0) const;
  Scalar computeA(const Scalar& R, const Scalar& alpha) const;
  Scalar computeH(const Scalar& R, const Scalar& alpha) const;
  Scalar computeApproxA(const Scalar& alpha, const Scalar& d0) const;
  Scalar computeApproxdEdd(const Scalar& alpha, const Scalar& d0) const;

  Scalar computeRPlanar(const Scalar& alpha, const Scalar& d0) const;
  Scalar computeAPlanar(const Scalar& R, const Scalar& alpha) const;
  Scalar computeHPlanar(const Scalar& R, const Scalar& alpha) const;
  Scalar computeApproxAPlanar(const Scalar& alpha, const Scalar& d0) const;
  Scalar computeApproxdEddPlanar(const Scalar& alpha, const Scalar& d0) const;

  Scalar computedEdd(const Scalar& R, const Scalar& alpha) const;
  Scalar computedEddPlanar(const Scalar& R, const Scalar& alpha) const;

  Scalar interpolate_table(const Scalar& A, const Scalar& d0, const MatXx& mat,
                           const Scalar& dmin) const;
  Scalar interpolate_table_planar(const Scalar& A, const Scalar& d0,
                                  const MatXx& mat, const Scalar& dmin) const;
  Scalar interpolate_table_grad(const Scalar& A, const Scalar& d0,
                                const MatXx& mat, const Scalar& dmin) const;
  Scalar interpolate_table_grad_planar(const Scalar& A, const Scalar& d0,
                                       const MatXx& mat,
                                       const Scalar& dmin) const;

  Scalar computedEddAreaDist(const Scalar& A_target, const Scalar& d0) const;
  Scalar computedEddAreaDistPlanar(const Scalar& A_target,
                                   const Scalar& d0) const;

  void print_energy_data(std::ostream& oss, bool first_time) const;
  void print_table(std::ostream& oss, const MatXx& mat,
                   const Scalar& dmin) const;
  void print_dEdd_table(std::ostream& oss) const;

  CohesionTable();

  void setParameter(const Scalar& sigma, const Scalar& theta,
                    const Scalar& radii, const Scalar& max_d0,
                    const int disc = 256);

  void construct_alpha_table();
  void construct_planar_alpha_table();

  Scalar getRadii() const;

  Scalar getDMin() const;

  Scalar getDMinPlanar() const;

  Scalar interpolate_dEdd(const Scalar& A, const Scalar& d0) const;
  Scalar interpolate_d2Edd2(const Scalar& A, const Scalar& d0) const;
  Scalar interpolate_alpha(const Scalar& A, const Scalar& d0) const;

  Scalar interpolate_dEdd_planar(const Scalar& A, const Scalar& d0) const;
  Scalar interpolate_d2Edd2_planar(const Scalar& A, const Scalar& d0) const;
  Scalar interpolate_alpha_planar(const Scalar& A, const Scalar& d0) const;
};
};  // namespace strandsim

#endif
