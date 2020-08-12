/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * So-bogus is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * So-bogus is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with So-bogus.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BOGUS_FISCHER_BURMEISTER_IMPL_HPP
#define BOGUS_FISCHER_BURMEISTER_IMPL_HPP

#include <iostream>

#include "FischerBurmeister.hpp"
namespace bogus {

template <DenseIndexType Dimension, typename Scalar>
Scalar SimpleFBBaseFunction<Dimension, Scalar>::compute(const Scalar& x,
                                                        const Scalar& y) {
  return x + y - std::sqrt(x * x + y * y);
}

template <DenseIndexType Dimension, typename Scalar>
Scalar SimpleFBBaseFunction<Dimension, Scalar>::computeJacobian(
    const Scalar& x, const Scalar& y, Scalar& dFb_dx, Scalar& dFb_dy) {
  if (x == 0. && y == 0.) {
    dFb_dx = 0.;
    dFb_dy = 0.;
  } else {
    dFb_dx = 1 - x / std::sqrt(x * x + y * y);
    dFb_dy = 1 - y / std::sqrt(x * x + y * y);
  }

  return x + y - std::sqrt(x * x + y * y);
}

template <DenseIndexType Dimension, typename Scalar>
Scalar SimpleFBBaseFunction<Dimension, Scalar>::computeHessian(
    const Scalar& x, const Scalar& y, Scalar& d2Fb_dx2, Scalar& d2Fb_dxy,
    Scalar& d2Fb_dy2) {
  if (x == 0. && y == 0.) {
    d2Fb_dx2 = 0.;
    d2Fb_dy2 = 0.;
    d2Fb_dxy = 0.;
  } else {
    const Scalar t0 = 1. / pow(x * x + y * y, 1.5);
    d2Fb_dx2 = -y * y * t0;
    d2Fb_dy2 = -x * x * t0;
    d2Fb_dxy = x * y * t0;
  }

  return 0.;
}

template <DenseIndexType Dimension, typename Scalar>
template <bool JacobianAsWell>
Scalar SimpleFBBaseFunction<Dimension, Scalar>::compute(const Scalar& x,
                                                        const Scalar& y,
                                                        Scalar& dFb_dx,
                                                        Scalar& dFb_dy) {
  if (JacobianAsWell) {
    if (x == 0. && y == 0.) {
      dFb_dx = 0.;
      dFb_dy = 0.;
    } else {
      dFb_dx = 1 - x / std::sqrt(x * x + y * y);
      dFb_dy = 1 - y / std::sqrt(x * x + y * y);
    }
  }

  return x + y - std::sqrt(x * x + y * y);
}

template <DenseIndexType Dimension, typename Scalar>
void FBBaseFunction<Dimension, Scalar>::compute(const Scalar mu,
                                                const Vector& x,
                                                const Vector& y, Vector& fb) {
  static Matrix a, b;  // unused

  if (NumTraits<Scalar>::isZero(mu)) {
    Traits::np(fb) = x[0] + y[0] - std::sqrt(x[0] * x[0] + y[0] * y[0]);
    Traits::tp(fb) = Traits::tp(x);
  } else {
    Vector xh = x, yh = y;
    Traits::np(xh) *= mu;
    Traits::tp(yh) *= mu;
    compute<false>(xh, yh, fb, a, b);
  }
}

template <DenseIndexType Dimension, typename Scalar>
void FBBaseFunction<Dimension, Scalar>::computeJacobian(
    const Scalar mu, const Vector& x, const Vector& y, Vector& fb,
    Matrix& dFb_dx, Matrix& dFb_dy) {
  if (NumTraits<Scalar>::isZero(mu)) {
    const Scalar z = std::sqrt(x[0] * x[0] + y[0] * y[0]);
    Traits::np(fb) = x[0] + y[0] - z;
    Traits::tp(fb) = Traits::tp(x);

    dFb_dx.setIdentity();
    dFb_dy.setZero();

    if (!NumTraits<Scalar>::isZero(z)) {
      dFb_dx(0, 0) = 1. - x[0] / z;
      dFb_dy(0, 0) = 1. - y[0] / z;
    }
  } else {
    Vector xh = x, yh = y;
    Traits::np(xh) *= mu;
    Traits::tp(yh) *= mu;
    compute<true>(xh, yh, fb, dFb_dx, dFb_dy);
    Traits::nc(dFb_dx) *= mu;
    Traits::tc(dFb_dy) *= mu;
  }
}

template <DenseIndexType Dimension, typename Scalar>
template <bool JacobianAsWell>
void FBBaseFunction<Dimension, Scalar>::compute(const Vector& x,
                                                const Vector& y, Vector& fb,
                                                Matrix& dFb_dx,
                                                Matrix& dFb_dy) {
  const unsigned d = x.rows();

  // see [Daviet et al 2011], Appendix A.1

  Vector z2(d);

  z2[0] = x.squaredNorm() + y.squaredNorm();
  Traits::tp(z2) = x[0] * Traits::tp(x) + y[0] * Traits::tp(y);
  const Scalar nz2t = Traits::tp(z2).norm();

  Vector omega1(d), omega2(d);
  omega1[0] = omega2[0] = .5;

  if (NumTraits<Scalar>::isZero(nz2t)) {
    Traits::tp(omega1).setZero();
    omega1[1] = -.5;
    Traits::tp(omega2).setZero();
    omega2[1] = .5;
  } else {
    Traits::tp(omega1) = -.5 * Traits::tp(z2) / nz2t;
    Traits::tp(omega2) = .5 * Traits::tp(z2) / nz2t;
  }

  const Scalar rlambda1 = std::sqrt(std::max((Scalar)0, z2[0] - 2 * nz2t));
  const Scalar rlambda2 = std::sqrt(z2[0] + 2 * nz2t);

  const Vector z = rlambda1 * omega1 + rlambda2 * omega2;
  fb = x + y - z;

  if (JacobianAsWell) {
    const Matrix Id(Matrix::Identity(d, d));

    if (NumTraits<Scalar>::isZero(rlambda2)) {
      // x = y = 0
      dFb_dx.setZero(d, d);
      dFb_dy.setZero(d, d);
    } else {
      if (NumTraits<Scalar>::isZero(rlambda1)) {
        const Scalar izn = 1. / (x[0] * x[0] + y[0] * y[0]);
        dFb_dx = (1. - x[0] * izn) * Id;
        dFb_dy = (1. - y[0] * izn) * Id;
      } else {
        const Scalar det = rlambda1 * rlambda2;

        Matrix L, invLz(d, d);

        invLz(0, 0) = z[0];
        Traits::ntb(invLz) = -Traits::tp(z).transpose();
        Traits::tnb(invLz) = -Traits::tp(z);
        Traits::ttb(invLz) = (det * Traits::TgMatrix::Identity(d - 1, d - 1) +
                              Traits::tp(z) * Traits::tp(z).transpose()) /
                             z[0];
        invLz /= det;

        L = x[0] * Id;
        Traits::ntb(L) = Traits::tp(x).transpose();
        Traits::tnb(L) = Traits::tp(x);
        dFb_dx.setIdentity(d, d);
        dFb_dx.noalias() -= invLz * L;

        L = y[0] * Id;
        Traits::ntb(L) = Traits::tp(y).transpose();
        Traits::tnb(L) = Traits::tp(y);
        dFb_dy.setIdentity(d, d);
        dFb_dy.noalias() -= invLz * L;
      }
    }
  }
}

template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV>
void FischerBurmeister<Dimension, Scalar, DeSaxceCOV>::compute(const Scalar mu,
                                                               const Vector& x,
                                                               const Vector& y,
                                                               Vector& fb) {
  if (DeSaxceCOV) {
    Vector yt(y);
    Traits::np(yt) += mu * Traits::tp(y).norm();
    BaseFunction::compute(mu, x, yt, fb);
  } else {
    BaseFunction::compute(mu, x, y, fb);
  }
}

template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV>
void FischerBurmeister<Dimension, Scalar, DeSaxceCOV>::compute(
    const Vector& x, Vector& fb) const {
  const Vector y = m_A * x + m_b;
  const Vector xs = m_scaling * x;
  compute(m_mu, xs, y, fb);
}

template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV>
void FischerBurmeister<Dimension, Scalar, DeSaxceCOV>::computeJacobian(
    const Vector& x, Vector& fb, Matrix& dFb_dx) const {
  const Vector xs = m_scaling * x;
  Vector y(m_A * x + m_b);
  Scalar s = 0.;
  //  std::cout << m_A << std::endl ;
  //  std::cout << m_b.transpose() << std::endl ;

  if (DeSaxceCOV) {
    s = Traits::tp(y).norm();
    Traits::np(y) += m_mu * s;
  }
  //  std::cout << y.transpose() << std::endl ;

  Matrix dFb_dy;
  BaseFunction::computeJacobian(m_mu, xs, y, fb, dFb_dx, dFb_dy);
  //  std::cout << dFb_dx.transpose() << std::endl ;
  //  std::cout << dFb_dy.transpose() << std::endl ;

  if (DeSaxceCOV && !NumTraits<Scalar>::isZero(s)) {
    Traits::tc(dFb_dy).noalias() +=
        Traits::nc(dFb_dy) * (m_mu / s) * Traits::tp(y).transpose();
  }

  dFb_dx *= m_scaling;
  dFb_dx.noalias() += dFb_dy * m_A;
}

template <DenseIndexType Dimension, typename Scalar>
void HerschelBulkleyFischerBurmeister<Dimension, Scalar>::compute(
    const Vector& x, Vector& fb) const {
  Vector y = m_A * x + m_b;
  compute(m_yield, m_eta, m_n, m_scaling, x, y, m_A, fb);
}

template <DenseIndexType Dimension, typename Scalar>
void HerschelBulkleyFischerBurmeister<Dimension, Scalar>::computeJacobian(
    const Vector& x, Vector& fb, Matrix& dFb_dx) const {
  Vector y = m_A * x + m_b;

  compute(m_yield, m_eta, m_n, m_scaling, x, y, m_A, fb);
  computeJacobian(m_yield, m_eta, m_n, m_scaling, x, y, m_A, dFb_dx);
}

template <DenseIndexType Dimension, typename Scalar>
void HerschelBulkleyFischerBurmeister<Dimension, Scalar>::d2xdr2(
    const Scalar eta, const Scalar n, const Scalar scaling, const Vector& x,
    const Vector& y, const Matrix& A, Matrix& Matd2xdr2) {
  Mat2x W = A.template block<Dimension - 1, Dimension - 1>(1, 1);
  Vec2x v = Traits::tp(y);
  Vec2x r = Traits::tp(x);
  const Scalar nv = v.norm();
  const Scalar nr = r.norm();

  Matd2xdr2.setZero();

  if (nv < 1e-20 || nr < 1e-20) {
    return;
  }

  Vec2x vhat = v / nv;
  Vec2x rhat = r / nr;

  Traits::ttb(Matd2xdr2) =
      scaling *
      (eta * n *
           ((n - 1.) * pow(nv, n - 2.) * W.transpose() * vhat *
                vhat.transpose() * W +
            pow(nv, n - 1.) * W.transpose() *
                (Mat2x::Identity() - vhat * vhat.transpose()) / nv * W) -
       (Mat2x::Identity() - rhat * rhat.transpose()) / nr);
}

template <DenseIndexType Dimension, typename Scalar>
void HerschelBulkleyFischerBurmeister<Dimension, Scalar>::d2ydr2(
    const Vector& y, const Matrix& A, Matrix& Matd2ydr2) {
  Mat2x W = A.template block<Dimension - 1, Dimension - 1>(1, 1);
  Vec2x v = Traits::tp(y);

  const Scalar nv = v.norm();

  Matd2ydr2.setZero();

  if (nv < 1e-20) {
    return;
  }

  Vec2x vhat = v / nv;
  Traits::ttb(Matd2ydr2) =
      W.transpose() * (Mat2x::Identity() - vhat * vhat.transpose()) / nv * W;
}

template <DenseIndexType Dimension, typename Scalar>
void HerschelBulkleyFischerBurmeister<Dimension, Scalar>::computeJacobian(
    const Scalar yield, const Scalar eta, const Scalar n, const Scalar scaling,
    const Vector& x, const Vector& y, const Matrix& A, Matrix& fb) {
  Vector xdr, ydr;

  dxdr(eta, n, scaling, x, y, A, xdr);
  dydr(y, A, ydr);

  Matrix x2dr2, y2dr2;
  d2xdr2(eta, n, scaling, x, y, A, x2dr2);
  d2ydr2(y, A, y2dr2);

  const Scalar ny = Traits::tp(y).norm();
  const Scalar nx = Traits::tp(x).norm();

  const Scalar rTA = yield + eta * pow(ny, n) - nx;

  Scalar d2Fbdx2(0), d2Fbdxy(0), d2Fbdy2(0);
  BaseFunction::computeHessian(rTA * scaling, ny, d2Fbdx2, d2Fbdxy, d2Fbdy2);

  Scalar dfdx(0), dfdy(0);
  BaseFunction::computeJacobian(rTA * scaling, ny, dfdx, dfdy);

  fb = d2Fbdx2 * xdr * xdr.transpose() +
       d2Fbdxy * (xdr * ydr.transpose() + ydr * xdr.transpose()) +
       d2Fbdy2 * ydr * ydr.transpose() + dfdx * x2dr2 + dfdy * y2dr2;
  fb(0, 0) = 1.0;
}

template <DenseIndexType Dimension, typename Scalar>
void HerschelBulkleyFischerBurmeister<Dimension, Scalar>::compute(
    const Scalar yield, const Scalar eta, const Scalar n, const Scalar scaling,
    const Vector& x, const Vector& y, const Matrix& A, Vector& fb) {
  const Scalar ny = Traits::tp(y).norm();
  const Scalar nx = Traits::tp(x).norm();

  const Scalar rTA = yield + eta * pow(ny, n) - nx;
  const Scalar f = BaseFunction::compute(scaling * rTA, ny);

  Scalar dfdx(0), dfdy(0);
  BaseFunction::computeJacobian(rTA * scaling, ny, dfdx, dfdy);
  Vector xdr, ydr;

  dxdr(eta, n, scaling, x, y, A, xdr);
  dydr(y, A, ydr);

  Vector gfb = dfdx * xdr + dfdy * ydr;
  fb = f * gfb;
}

template <DenseIndexType Dimension, typename Scalar>
Scalar HerschelBulkleyFischerBurmeister<Dimension, Scalar>::func(
    const Scalar yield, const Scalar eta, const Scalar n, const Scalar scaling,
    const Vector& x, const Vector& y) {
  const Scalar ny = Traits::tp(y).norm();
  const Scalar nx = Traits::tp(x).norm();

  const Scalar rTA = yield + eta * pow(ny, n) - nx;
  return BaseFunction::compute(scaling * rTA, ny);
}

template <DenseIndexType Dimension, typename Scalar>
void HerschelBulkleyFischerBurmeister<Dimension, Scalar>::dxdr(
    const Scalar eta, const Scalar n, const Scalar scaling, const Vector& x,
    const Vector& y, const Matrix& A, Vector& fb) {
  const Scalar ny = Traits::tp(y).norm();
  const Scalar nx = Traits::tp(x).norm();
  fb.setZero();

  Mat2x W = A.template block<Dimension - 1, Dimension - 1>(1, 1);

  if (ny < 1e-20 || nx < 1e-20) {
    // shouldn't be here for nx actually
    return;
  } else {
    const Vector nnx = x / nx;
    const Vector nnv = y / ny;
    Traits::tp(fb) =
        (eta * n * pow(ny, n - 1.) * W.transpose() * Traits::tp(nnv) -
         Traits::tp(nnx)) *
        scaling;
  }
}

template <DenseIndexType Dimension, typename Scalar>
void HerschelBulkleyFischerBurmeister<Dimension, Scalar>::dydr(const Vector& y,
                                                               const Matrix& A,
                                                               Vector& fb) {
  const Scalar ny = Traits::tp(y).norm();
  fb.setZero();
  if (ny < 1e-20) {
    return;
  } else {
    Mat2x W = A.template block<Dimension - 1, Dimension - 1>(1, 1);

    Traits::tp(fb) = W.transpose() * Traits::tp(y) / ny;
  }
}

template <DenseIndexType Dimension, typename Scalar>
Scalar HerschelBulkleyFischerBurmeister<Dimension, Scalar>::yieldFunction(
    const Scalar yield, const Scalar eta, const Scalar n, const Vector& x,
    const Vector& y) {
  const Scalar ny = Traits::tp(y).norm();
  const Scalar nx = Traits::tp(x).norm();

  const Scalar rTA = yield + eta * pow(ny, n) - nx;
  return rTA;
}

template <DenseIndexType Dimension, typename Scalar>
Scalar HerschelBulkleyFischerBurmeister<Dimension, Scalar>::yieldFunction(
    const Vector& x) {
  Vector y = m_A * x + m_b;

  return yieldFunction(m_yield, m_eta, m_n, x, y);
}

template <DenseIndexType Dimension, typename Scalar>
Scalar HerschelBulkleyFischerBurmeister<Dimension, Scalar>::computeEnergy(
    const Scalar yield, const Scalar eta, const Scalar n, const Scalar scaling,
    const Vector& x, const Vector& y) {
  const Scalar ny = Traits::tp(y).norm();
  const Scalar nx = Traits::tp(x).norm();

  const Scalar rTA = yield + eta * pow(ny, n) - nx;
  const Scalar f = BaseFunction::compute(scaling * rTA, ny);

  return 0.5 * f * f;
}
}  // namespace bogus

#endif
