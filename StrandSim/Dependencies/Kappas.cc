/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Kappas.hh"

#include "../Core/ElasticStrandUtils.hh"
#include "../Utils/MathUtilities.hh"

//#define CURVATURE_USE_SINTHETA

#ifndef CURVATURE_USE_SINTHETA
#define CURVATURE_USE_TANTHETA
#endif

namespace strandsim {

void Kappas::compute() {
  m_value.resize(m_size);
  const Vec3xArray& curvatureBinormals = m_curvatureBinormals.get();
  const Vec3xArray& materialFrames1 = m_materialFrames1.get();
  const Vec3xArray& materialFrames2 = m_materialFrames2.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    const Vec3x& kb = curvatureBinormals[vtx];
    const Vec3x& m1e = materialFrames1[vtx - 1];
    const Vec3x& m2e = materialFrames2[vtx - 1];
    const Vec3x& m1f = materialFrames1[vtx];
    const Vec3x& m2f = materialFrames2[vtx];

    m_value[vtx] = Vec4x(kb.dot(m2e), -kb.dot(m1e), kb.dot(m2f), -kb.dot(m1f));
  }

  setDependentsDirty();
}

void GradKappas::compute() {
  m_value.resize(m_size);
  const std::vector<Scalar>& lengths = m_lengths.get();
  const Vec3xArray& tangents = m_tangents.get();
  const Vec3xArray& curvatureBinormals = m_curvatureBinormals.get();
  const Vec3xArray& materialFrames1 = m_materialFrames1.get();
  const Vec3xArray& materialFrames2 = m_materialFrames2.get();
  const Vec4xArray& kappas = m_kappas.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    GradKType& gradKappa = m_value[vtx];

    const Scalar norm_e = lengths[vtx - 1];
    const Scalar norm_f = lengths[vtx];

    const Vec3x& te = tangents[vtx - 1];
    const Vec3x& tf = tangents[vtx];

    const Vec3x& m1e = materialFrames1[vtx - 1];
    const Vec3x& m2e = materialFrames2[vtx - 1];
    const Vec3x& m1f = materialFrames1[vtx];
    const Vec3x& m2f = materialFrames2[vtx];
#ifdef CURVATURE_USE_SINTHETA
    Scalar chi = (te + tf).norm();

    if (chi < 1e-12) {
      WarningStream(g_log, "")
          << "GradKappas::compute(): "
          << " chi = " << chi << " te = " << te << " tf = " << tf;
      chi = 1e-12;
    }

    const Vec3x& tilde_t = (te + tf) / chi;

    const Vec4x& kappa = kappas[vtx];

    const Vec3x& bar_te = tilde_t + tilde_t.dot(tf) * te;
    const Vec3x& bar_tf = tilde_t + tilde_t.dot(te) * tf;

    const Vec3x& Dkappa0eDe =
        1.0 / (norm_e * chi) * (-kappa[0] * bar_te + tf.cross(2.0 * m2e));
    const Vec3x& Dkappa0eDf =
        1.0 / (norm_f * chi) * (-kappa[0] * bar_tf - te.cross(2.0 * m2e));
    const Vec3x& Dkappa1eDe =
        1.0 / (norm_e * chi) * (-kappa[1] * bar_te - tf.cross(2.0 * m1e));
    const Vec3x& Dkappa1eDf =
        1.0 / (norm_f * chi) * (-kappa[1] * bar_tf + te.cross(2.0 * m1e));

    const Vec3x& Dkappa0fDe =
        1.0 / (norm_e * chi) * (-kappa[2] * bar_te + tf.cross(2.0 * m2f));
    const Vec3x& Dkappa0fDf =
        1.0 / (norm_f * chi) * (-kappa[2] * bar_tf - te.cross(2.0 * m2f));
    const Vec3x& Dkappa1fDe =
        1.0 / (norm_e * chi) * (-kappa[3] * bar_te - tf.cross(2.0 * m1f));
    const Vec3x& Dkappa1fDf =
        1.0 / (norm_f * chi) * (-kappa[3] * bar_tf + te.cross(2.0 * m1f));
#else
    Scalar chi = 1.0 + te.dot(tf);

    //    assert( chi>0 );
    if (chi < 1e-12) {
      WarningStream(g_log, "")
          << "GradKappas::compute(): "
          << " chi = " << chi << " te = " << te << " tf = " << tf;
      chi = 1e-12;
    }

    const Vec3x& tilde_t = (te + tf) / chi;
    const Vec3x& tilde_d1e = (2.0 * m1e) / chi;
    const Vec3x& tilde_d1f = (2.0 * m1f) / chi;
    const Vec3x& tilde_d2e = (2.0 * m2e) / chi;
    const Vec3x& tilde_d2f = (2.0 * m2f) / chi;

    const Vec4x& kappa = kappas[vtx];

    const Vec3x& Dkappa0eDe =
        1.0 / norm_e * (-kappa[0] * tilde_t + tf.cross(tilde_d2e));
    const Vec3x& Dkappa0eDf =
        1.0 / norm_f * (-kappa[0] * tilde_t - te.cross(tilde_d2e));
    const Vec3x& Dkappa1eDe =
        1.0 / norm_e * (-kappa[1] * tilde_t - tf.cross(tilde_d1e));
    const Vec3x& Dkappa1eDf =
        1.0 / norm_f * (-kappa[1] * tilde_t + te.cross(tilde_d1e));

    const Vec3x& Dkappa0fDe =
        1.0 / norm_e * (-kappa[2] * tilde_t + tf.cross(tilde_d2f));
    const Vec3x& Dkappa0fDf =
        1.0 / norm_f * (-kappa[2] * tilde_t - te.cross(tilde_d2f));
    const Vec3x& Dkappa1fDe =
        1.0 / norm_e * (-kappa[3] * tilde_t - tf.cross(tilde_d1f));
    const Vec3x& Dkappa1fDf =
        1.0 / norm_f * (-kappa[3] * tilde_t + te.cross(tilde_d1f));
#endif
    gradKappa.block<3, 1>(0, 0) = -Dkappa0eDe;
    gradKappa.block<3, 1>(4, 0) = Dkappa0eDe - Dkappa0eDf;
    gradKappa.block<3, 1>(8, 0) = Dkappa0eDf;
    gradKappa.block<3, 1>(0, 1) = -Dkappa1eDe;
    gradKappa.block<3, 1>(4, 1) = Dkappa1eDe - Dkappa1eDf;
    gradKappa.block<3, 1>(8, 1) = Dkappa1eDf;

    gradKappa.block<3, 1>(0, 2) = -Dkappa0fDe;
    gradKappa.block<3, 1>(4, 2) = Dkappa0fDe - Dkappa0fDf;
    gradKappa.block<3, 1>(8, 2) = Dkappa0fDf;
    gradKappa.block<3, 1>(0, 3) = -Dkappa1fDe;
    gradKappa.block<3, 1>(4, 3) = Dkappa1fDe - Dkappa1fDf;
    gradKappa.block<3, 1>(8, 3) = Dkappa1fDf;

    const Vec3x& kb = curvatureBinormals[vtx];

    gradKappa(3, 0) = -kb.dot(m1e);
    gradKappa(7, 0) = 0.0;
    gradKappa(3, 1) = -kb.dot(m2e);
    gradKappa(7, 1) = 0.0;

    gradKappa(3, 2) = 0.0;
    gradKappa(7, 2) = -kb.dot(m1f);
    gradKappa(3, 3) = 0.0;
    gradKappa(7, 3) = -kb.dot(m2f);
  }

  setDependentsDirty();
}

void HessKappas::compute() {
  m_value.resize(m_size);
  const std::vector<Scalar>& lengths = m_lengths.get();
  const Vec3xArray& tangents = m_tangents.get();
  const Vec3xArray& curvatureBinormals = m_curvatureBinormals.get();
  const Vec3xArray& materialFrames1 = m_materialFrames1.get();
  const Vec3xArray& materialFrames2 = m_materialFrames2.get();
  const Vec4xArray& kappas = m_kappas.get();

  for (IndexType vtx = m_firstValidIndex; vtx < m_curvatureBinormals.size();
       ++vtx) {
    Mat11x& DDkappa0 = m_value[vtx * 4 + 0];
    Mat11x& DDkappa1 = m_value[vtx * 4 + 1];
    Mat11x& DDkappa2 = m_value[vtx * 4 + 2];
    Mat11x& DDkappa3 = m_value[vtx * 4 + 3];

    DDkappa0.setZero();
    DDkappa1.setZero();
    DDkappa2.setZero();
    DDkappa3.setZero();

    const Scalar norm_e = lengths[vtx - 1];
    const Scalar norm_f = lengths[vtx];
    const Scalar norm2_e = square(
        norm_e);  // That's bloody stupid, taking the square of a square root.
    const Scalar norm2_f = square(norm_f);

    const Vec3x& te = tangents[vtx - 1];
    const Vec3x& tf = tangents[vtx];

    const Vec3x& m1e = materialFrames1[vtx - 1];
    const Vec3x& m2e = materialFrames2[vtx - 1];
    const Vec3x& m1f = materialFrames1[vtx];
    const Vec3x& m2f = materialFrames2[vtx];
#ifdef CURVATURE_USE_SINTHETA
    Scalar chi = (te + tf).norm();

    if (chi < 1e-12) {
      WarningStream(g_log, "")
          << "HessKappas::compute(): "
          << " chi = " << chi << " te = " << te << " tf = " << tf;
      chi = 1e-12;
    }

    const Vec3x& tilde_t = (te + tf) / chi;

    const Vec4x& kappa = kappas[vtx];

    const Vec3x& bar_te = tilde_t + tilde_t.dot(tf) * te;
    const Vec3x& bar_tf = tilde_t + tilde_t.dot(te) * tf;

    const Vec3x& Dkappa0eDe =
        1.0 / (norm_e * chi) * (-kappa[0] * bar_te + tf.cross(2.0 * m2e));
    const Vec3x& Dkappa0eDf =
        1.0 / (norm_f * chi) * (-kappa[0] * bar_tf - te.cross(2.0 * m2e));
    const Vec3x& Dkappa1eDe =
        1.0 / (norm_e * chi) * (-kappa[1] * bar_te - tf.cross(2.0 * m1e));
    const Vec3x& Dkappa1eDf =
        1.0 / (norm_f * chi) * (-kappa[1] * bar_tf + te.cross(2.0 * m1e));

    const Vec3x& Dkappa0fDe =
        1.0 / (norm_e * chi) * (-kappa[2] * bar_te + tf.cross(2.0 * m2f));
    const Vec3x& Dkappa0fDf =
        1.0 / (norm_f * chi) * (-kappa[2] * bar_tf - te.cross(2.0 * m2f));
    const Vec3x& Dkappa1fDe =
        1.0 / (norm_e * chi) * (-kappa[3] * bar_te - tf.cross(2.0 * m1f));
    const Vec3x& Dkappa1fDf =
        1.0 / (norm_f * chi) * (-kappa[3] * bar_tf + te.cross(2.0 * m1f));

    const Vec3x& kb = curvatureBinormals[vtx];

    const Mat3x& Id = Mat3x::Identity();

    const Mat3x& DteDe = 1.0 / norm_e * (Id - outerProd<3>(te, te));
    const Mat3x& DtfDf = 1.0 / norm_f * (Id - outerProd<3>(tf, tf));

    const Vec3x& DchiDe = DteDe * tilde_t;
    const Vec3x& DchiDf = DtfDf * tilde_t;

    const Mat3x& DttDe = 1.0 / (norm_e * chi) *
                         (Id - outerProd<3>(tilde_t, tilde_t)) *
                         (Id - outerProd<3>(te, te));
    const Mat3x& DttDf = 1.0 / (norm_f * chi) *
                         (Id - outerProd<3>(tilde_t, tilde_t)) *
                         (Id - outerProd<3>(tf, tf));

    const Mat3x& DbteDe =
        (Id + outerProd<3>(te, tf)) * DttDe + tilde_t.dot(tf) * DteDe;
    const Mat3x& DbteDf =
        (Id + outerProd<3>(te, tf)) * DtfDf + outerProd<3>(te, tilde_t) * DtfDf;

    const Mat3x& DbtfDf =
        (Id + outerProd<3>(tf, te)) * DttDf + tilde_t.dot(te) * DtfDf;

    // 1st Hessian
    const Mat3x& D2kappa0De2 =
        -1.0 / (norm_e * chi) *
        symPart<3>(outerProd<3>(bar_te + chi * te, Dkappa0eDe) +
                   outerProd<3>(Dkappa0eDe, DchiDe * norm_e) +
                   kappa[0] * DbteDe);

    const Mat3x& D2kappa0Df2 =
        -1.0 / (norm_f * chi) *
        symPart<3>(outerProd<3>(bar_tf + chi * tf, Dkappa0eDf) +
                   outerProd<3>(Dkappa0eDf, DchiDf * norm_f) +
                   kappa[0] * DbtfDf);

    const Mat3x& D2kappa0DeDf =
        -1.0 / (norm_e * chi) *
        (outerProd<3>(bar_te, Dkappa0eDf) +
         outerProd(Dkappa0eDe, DchiDf) * norm_e + kappa[0] * DbteDf +
         2.0 * crossMat(m2e) * DtfDf);

    const Mat3x& D2kappa0DfDe = D2kappa0DeDf.transpose();

    const Scalar D2kappa0Dthetae2 = -kb.dot(m2e);
    const Scalar D2kappa0Dthetaf2 = 0.0;
    const Vec3x& D2kappa0DeDthetae =
        1.0 / (chi * norm_e) * (kb.dot(m1e) * bar_te - tf.cross(2.0 * m1e));
    const Vec3x& D2kappa0DeDthetaf = Vec3x::Zero();
    const Vec3x& D2kappa0DfDthetae =
        1.0 / (chi * norm_f) * (kb.dot(m1e) * bar_tf + te.cross(2.0 * m1e));
    const Vec3x& D2kappa0DfDthetaf = Vec3x::Zero();

    // 2nd Hessian
    const Mat3x& D2kappa1De2 =
        -1.0 / (norm_e * chi) *
        symPart<3>(outerProd<3>(bar_te + chi * te, Dkappa1eDe) +
                   outerProd<3>(Dkappa1eDe, DchiDe * norm_e) +
                   kappa[1] * DbteDe);

    const Mat3x& D2kappa1Df2 =
        -1.0 / (norm_f * chi) *
        symPart<3>(outerProd<3>(bar_tf + chi * tf, Dkappa1eDf) +
                   outerProd<3>(Dkappa1eDf, DchiDf * norm_f) +
                   kappa[1] * DbtfDf);

    const Mat3x& D2kappa1DeDf =
        -1.0 / (norm_e * chi) *
        (outerProd<3>(bar_te, Dkappa1eDf) +
         outerProd(Dkappa1eDe, DchiDf) * norm_e + kappa[1] * DbteDf -
         2.0 * crossMat(m1e) * DtfDf);

    const Mat3x& D2kappa1DfDe = D2kappa1DeDf.transpose();

    const Scalar D2kappa1Dthetae2 = kb.dot(m1e);
    const Scalar D2kappa1Dthetaf2 = 0.0;
    const Vec3x& D2kappa1DeDthetae =
        1.0 / (chi * norm_e) * (kb.dot(m2e) * bar_te - tf.cross(2.0 * m2e));
    const Vec3x& D2kappa1DeDthetaf = Vec3x::Zero();
    const Vec3x& D2kappa1DfDthetae =
        1.0 / (chi * norm_f) * (kb.dot(m2e) * bar_tf + te.cross(2.0 * m2e));
    const Vec3x& D2kappa1DfDthetaf = Vec3x::Zero();

    // 3rd Hessian
    const Mat3x& D2kappa2De2 =
        -1.0 / (norm_e * chi) *
        symPart<3>(outerProd<3>(bar_te + chi * te, Dkappa0fDe) +
                   outerProd<3>(Dkappa0fDe, DchiDe * norm_e) +
                   kappa[2] * DbteDe);

    const Mat3x& D2kappa2Df2 =
        -1.0 / (norm_f * chi) *
        symPart<3>(outerProd<3>(bar_tf + chi * tf, Dkappa0fDf) +
                   outerProd<3>(Dkappa0fDf, DchiDf * norm_f) +
                   kappa[2] * DbtfDf);

    const Mat3x& D2kappa2DeDf =
        -1.0 / (norm_e * chi) *
        (outerProd<3>(bar_te, Dkappa0fDf) +
         outerProd(Dkappa0fDe, DchiDf) * norm_e + kappa[2] * DbteDf +
         2.0 * crossMat(m2f) * DtfDf);

    const Mat3x& D2kappa2DfDe = D2kappa2DeDf.transpose();

    const Scalar D2kappa2Dthetae2 = 0.0;
    const Scalar D2kappa2Dthetaf2 = -kb.dot(m2f);
    const Vec3x& D2kappa2DeDthetae = Vec3x::Zero();
    const Vec3x& D2kappa2DeDthetaf =
        1.0 / (chi * norm_e) * (kb.dot(m1f) * bar_te - tf.cross(2.0 * m1f));
    const Vec3x& D2kappa2DfDthetae = Vec3x::Zero();
    const Vec3x& D2kappa2DfDthetaf =
        1.0 / (chi * norm_f) * (kb.dot(m1f) * bar_tf + te.cross(2.0 * m1f));

    // 4th Hessian
    const Mat3x& D2kappa3De2 =
        -1.0 / (norm_e * chi) *
        symPart<3>(outerProd<3>(bar_te + chi * te, Dkappa1fDe) +
                   outerProd<3>(Dkappa1fDe, DchiDe * norm_e) +
                   kappa[3] * DbteDe);

    const Mat3x& D2kappa3Df2 =
        -1.0 / (norm_f * chi) *
        symPart<3>(outerProd<3>(bar_tf + chi * tf, Dkappa1fDf) +
                   outerProd<3>(Dkappa1fDf, DchiDf * norm_f) +
                   kappa[3] * DbtfDf);

    const Mat3x& D2kappa3DeDf =
        -1.0 / (norm_e * chi) *
        (outerProd<3>(bar_te, Dkappa1fDf) +
         outerProd(Dkappa1fDe, DchiDf) * norm_e + kappa[3] * DbteDf -
         2.0 * crossMat(m1f) * DtfDf);

    const Mat3x& D2kappa3DfDe = D2kappa3DeDf.transpose();

    const Scalar D2kappa3Dthetae2 = 0.0;
    const Scalar D2kappa3Dthetaf2 = kb.dot(m1f);
    const Vec3x& D2kappa3DeDthetae = Vec3x::Zero();
    const Vec3x& D2kappa3DeDthetaf =
        1.0 / (norm_e * chi) * (kb.dot(m2f) * bar_te - tf.cross(2.0 * m2f));
    const Vec3x& D2kappa3DfDthetae = Vec3x::Zero();
    const Vec3x& D2kappa3DfDthetaf =
        1.0 / (norm_f * chi) * (kb.dot(m2f) * bar_tf + te.cross(2.0 * m2f));
#else
    Scalar chi = 1.0 + te.dot(tf);

    //    assert( chi>0 );
    if (chi < 1e-12) {
      WarningStream(g_log, "")
          << "HessKappas::compute(): "
          << " chi = " << chi << " te = " << te << " tf = " << tf;
      chi = 1e-12;
    }

    const Vec3x& tilde_t = (te + tf) / chi;
    const Vec3x& tilde_d1e = (2.0 * m1e) / chi;
    const Vec3x& tilde_d2e = (2.0 * m2e) / chi;
    const Vec3x& tilde_d1f = (2.0 * m1f) / chi;
    const Vec3x& tilde_d2f = (2.0 * m2f) / chi;

    const Vec4x& kappa = kappas[vtx];

    const Vec3x& Dkappa0eDe =
        1.0 / norm_e * (-kappa[0] * tilde_t + tf.cross(tilde_d2e));
    const Vec3x& Dkappa0eDf =
        1.0 / norm_f * (-kappa[0] * tilde_t - te.cross(tilde_d2e));
    const Vec3x& Dkappa1eDe =
        1.0 / norm_e * (-kappa[1] * tilde_t - tf.cross(tilde_d1e));
    const Vec3x& Dkappa1eDf =
        1.0 / norm_f * (-kappa[1] * tilde_t + te.cross(tilde_d1e));

    const Vec3x& Dkappa0fDe =
        1.0 / norm_e * (-kappa[2] * tilde_t + tf.cross(tilde_d2f));
    const Vec3x& Dkappa0fDf =
        1.0 / norm_f * (-kappa[2] * tilde_t - te.cross(tilde_d2f));
    const Vec3x& Dkappa1fDe =
        1.0 / norm_e * (-kappa[3] * tilde_t - tf.cross(tilde_d1f));
    const Vec3x& Dkappa1fDf =
        1.0 / norm_f * (-kappa[3] * tilde_t + te.cross(tilde_d1f));

    const Vec3x& kb = curvatureBinormals[vtx];

    const Mat3x& Id = Mat3x::Identity();

    const Vec3x& DchiDe = 1.0 / norm_e * (Id - outerProd<3>(te, te)) * tf;
    const Vec3x& DchiDf = 1.0 / norm_f * (Id - outerProd<3>(tf, tf)) * te;

    const Mat3x& DtfDf = 1.0 / norm_f * (Id - outerProd<3>(tf, tf));

    const Mat3x& DttDe =
        1.0 / (chi * norm_e) *
        ((Id - outerProd<3>(te, te)) -
         outerProd<3>(tilde_t, (Id - outerProd<3>(te, te)) * tf));

    const Mat3x& DttDf =
        1.0 / (chi * norm_f) *
        ((Id - outerProd<3>(tf, tf)) -
         outerProd<3>(tilde_t, (Id - outerProd<3>(tf, tf)) * te));

    // 1st Hessian
    const Mat3x& D2kappa0De2 =
        -1.0 / norm_e *
        symPart<3>(outerProd<3>(tilde_t + te, Dkappa0eDe) + kappa[0] * DttDe +
                   outerProd<3>(1.0 / chi * tf.cross(tilde_d2e), DchiDe));

    const Mat3x& D2kappa0Df2 =
        -1.0 / norm_f *
        symPart<3>(outerProd<3>(tilde_t + tf, Dkappa0eDf) + kappa[0] * DttDf +
                   outerProd<3>(1.0 / chi * te.cross(tilde_d2e), DchiDf));

    const Mat3x& D2kappa0DeDf =
        -1.0 / norm_e *
        (outerProd<3>(tilde_t, Dkappa0eDf) + kappa[0] * DttDf +
         outerProd<3>(1.0 / chi * tf.cross(tilde_d2e), DchiDf) +
         crossMat(tilde_d2e) * DtfDf);

    const Mat3x& D2kappa0DfDe = D2kappa0DeDf.transpose();

    const Scalar D2kappa0Dthetae2 = -kb.dot(m2e);
    const Scalar D2kappa0Dthetaf2 = 0.0;
    const Vec3x& D2kappa0DeDthetae =
        1.0 / norm_e * (kb.dot(m1e) * tilde_t - 2.0 / chi * tf.cross(m1e));
    const Vec3x& D2kappa0DeDthetaf = Vec3x::Zero();
    const Vec3x& D2kappa0DfDthetae =
        1.0 / norm_f * (kb.dot(m1e) * tilde_t + 2.0 / chi * te.cross(m1e));
    const Vec3x& D2kappa0DfDthetaf = Vec3x::Zero();

    // 2nd Hessian
    const Mat3x& D2kappa1De2 =
        -1.0 / norm_e *
        symPart<3>(outerProd<3>(tilde_t + te, Dkappa1eDe) + kappa[1] * DttDe +
                   outerProd<3>(1.0 / chi * tf.cross(tilde_d1e), DchiDe));

    const Mat3x& D2kappa1Df2 =
        -1.0 / norm_f *
        symPart<3>(outerProd<3>(tilde_t + tf, Dkappa1eDf) + kappa[1] * DttDf +
                   outerProd<3>(1.0 / chi * te.cross(tilde_d1e), DchiDf));

    const Mat3x& D2kappa1DeDf =
        -1.0 / norm_e *
        (outerProd<3>(tilde_t, Dkappa1eDf) + kappa[1] * DttDf -
         outerProd<3>(1.0 / chi * tf.cross(tilde_d1e), DchiDf) -
         crossMat(tilde_d1e) * DtfDf);

    const Mat3x& D2kappa1DfDe = D2kappa1DeDf.transpose();

    const Scalar D2kappa1Dthetae2 = kb.dot(m1e);
    const Scalar D2kappa1Dthetaf2 = 0.0;
    const Vec3x& D2kappa1DeDthetae =
        1.0 / norm_e * (kb.dot(m2e) * tilde_t - 2.0 / chi * tf.cross(m2e));
    const Vec3x& D2kappa1DeDthetaf = Vec3x::Zero();
    const Vec3x& D2kappa1DfDthetae =
        1.0 / norm_f * (kb.dot(m2e) * tilde_t + 2.0 / chi * te.cross(m2e));
    const Vec3x& D2kappa1DfDthetaf = Vec3x::Zero();

    // 3rd Hessian
    const Mat3x& D2kappa2De2 =
        -1.0 / norm_e *
        symPart<3>(outerProd<3>(tilde_t + te, Dkappa0fDe) + kappa[2] * DttDe +
                   outerProd<3>(1.0 / chi * tf.cross(tilde_d2f), DchiDe));

    const Mat3x& D2kappa2Df2 =
        -1.0 / norm_f *
        symPart<3>(outerProd<3>(tilde_t + tf, Dkappa0fDf) + kappa[2] * DttDf +
                   outerProd<3>(1.0 / chi * te.cross(tilde_d2f), DchiDf));

    const Mat3x& D2kappa2DeDf =
        -1.0 / norm_e *
        (outerProd<3>(tilde_t, Dkappa0fDf) + kappa[2] * DttDf +
         outerProd<3>(1.0 / chi * tf.cross(tilde_d2f), DchiDf) +
         crossMat(tilde_d2f) * DtfDf);

    const Mat3x& D2kappa2DfDe = D2kappa2DeDf.transpose();

    const Scalar D2kappa2Dthetae2 = 0.0;
    const Scalar D2kappa2Dthetaf2 = -kb.dot(m2f);
    const Vec3x& D2kappa2DeDthetae = Vec3x::Zero();
    const Vec3x& D2kappa2DeDthetaf =
        1.0 / norm_e * (kb.dot(m1f) * tilde_t - 2.0 / chi * tf.cross(m1f));
    const Vec3x& D2kappa2DfDthetae = Vec3x::Zero();
    const Vec3x& D2kappa2DfDthetaf =
        1.0 / norm_f * (kb.dot(m1f) * tilde_t + 2.0 / chi * te.cross(m1f));

    // 4th Hessian
    const Mat3x& D2kappa3De2 =
        -1.0 / norm_e *
        symPart<3>(outerProd<3>(tilde_t + te, Dkappa1fDe) + kappa[3] * DttDe +
                   outerProd<3>(1.0 / chi * tf.cross(tilde_d1f), DchiDe));

    const Mat3x& D2kappa3Df2 =
        -1.0 / norm_f *
        symPart<3>(outerProd<3>(tilde_t + tf, Dkappa1fDf) + kappa[3] * DttDf +
                   outerProd<3>(1.0 / chi * te.cross(tilde_d1f), DchiDf));

    const Mat3x& D2kappa3DeDf =
        -1.0 / norm_e *
        (outerProd<3>(tilde_t, Dkappa1fDf) + kappa[3] * DttDf -
         outerProd<3>(1.0 / chi * tf.cross(tilde_d1f), DchiDf) -
         crossMat(tilde_d1f) * DtfDf);

    const Mat3x& D2kappa3DfDe = D2kappa3DeDf.transpose();

    const Scalar D2kappa3Dthetae2 = 0.0;
    const Scalar D2kappa3Dthetaf2 = kb.dot(m1f);
    const Vec3x& D2kappa3DeDthetae = Vec3x::Zero();
    const Vec3x& D2kappa3DeDthetaf =
        1.0 / norm_e * (kb.dot(m2f) * tilde_t - 2.0 / chi * tf.cross(m2f));
    const Vec3x& D2kappa3DfDthetae = Vec3x::Zero();
    const Vec3x& D2kappa3DfDthetaf =
        1.0 / norm_f * (kb.dot(m2f) * tilde_t + 2.0 / chi * te.cross(m2f));
#endif
    DDkappa0.block<3, 3>(0, 0) = D2kappa0De2;
    DDkappa0.block<3, 3>(0, 4) = -D2kappa0De2 + D2kappa0DeDf;
    DDkappa0.block<3, 3>(4, 0) = -D2kappa0De2 + D2kappa0DfDe;
    DDkappa0.block<3, 3>(4, 4) =
        D2kappa0De2 - (D2kappa0DeDf + D2kappa0DfDe) + D2kappa0Df2;
    DDkappa0.block<3, 3>(0, 8) = -D2kappa0DeDf;
    DDkappa0.block<3, 3>(8, 0) = -D2kappa0DfDe;
    DDkappa0.block<3, 3>(4, 8) = D2kappa0DeDf - D2kappa0Df2;
    DDkappa0.block<3, 3>(8, 4) = D2kappa0DfDe - D2kappa0Df2;
    DDkappa0.block<3, 3>(8, 8) = D2kappa0Df2;
    DDkappa0(3, 3) = D2kappa0Dthetae2;
    DDkappa0(7, 7) = D2kappa0Dthetaf2;
    DDkappa0(3, 7) = DDkappa0(7, 3) = 0.;
    DDkappa0.block<3, 1>(0, 3) = -D2kappa0DeDthetae;
    DDkappa0.block<1, 3>(3, 0) = DDkappa0.block<3, 1>(0, 3).transpose();
    DDkappa0.block<3, 1>(4, 3) = D2kappa0DeDthetae - D2kappa0DfDthetae;
    DDkappa0.block<1, 3>(3, 4) = DDkappa0.block<3, 1>(4, 3).transpose();
    DDkappa0.block<3, 1>(8, 3) = D2kappa0DfDthetae;
    DDkappa0.block<1, 3>(3, 8) = DDkappa0.block<3, 1>(8, 3).transpose();
    DDkappa0.block<3, 1>(0, 7) = -D2kappa0DeDthetaf;
    DDkappa0.block<1, 3>(7, 0) = DDkappa0.block<3, 1>(0, 7).transpose();
    DDkappa0.block<3, 1>(4, 7) = D2kappa0DeDthetaf - D2kappa0DfDthetaf;
    DDkappa0.block<1, 3>(7, 4) = DDkappa0.block<3, 1>(4, 7).transpose();
    DDkappa0.block<3, 1>(8, 7) = D2kappa0DfDthetaf;
    DDkappa0.block<1, 3>(7, 8) = DDkappa0.block<3, 1>(8, 7).transpose();

    assert(isSymmetric(DDkappa0));

    DDkappa1.block<3, 3>(0, 0) = D2kappa1De2;
    DDkappa1.block<3, 3>(0, 4) = -D2kappa1De2 + D2kappa1DeDf;
    DDkappa1.block<3, 3>(4, 0) = -D2kappa1De2 + D2kappa1DfDe;
    DDkappa1.block<3, 3>(4, 4) =
        D2kappa1De2 - (D2kappa1DeDf + D2kappa1DfDe) + D2kappa1Df2;
    DDkappa1.block<3, 3>(0, 8) = -D2kappa1DeDf;
    DDkappa1.block<3, 3>(8, 0) = -D2kappa1DfDe;
    DDkappa1.block<3, 3>(4, 8) = D2kappa1DeDf - D2kappa1Df2;
    DDkappa1.block<3, 3>(8, 4) = D2kappa1DfDe - D2kappa1Df2;
    DDkappa1.block<3, 3>(8, 8) = D2kappa1Df2;
    DDkappa1(3, 3) = D2kappa1Dthetae2;
    DDkappa1(7, 7) = D2kappa1Dthetaf2;
    DDkappa1(3, 7) = DDkappa1(7, 3) = 0.;
    DDkappa1.block<3, 1>(0, 3) = -D2kappa1DeDthetae;
    DDkappa1.block<1, 3>(3, 0) = DDkappa1.block<3, 1>(0, 3).transpose();
    DDkappa1.block<3, 1>(4, 3) = D2kappa1DeDthetae - D2kappa1DfDthetae;
    DDkappa1.block<1, 3>(3, 4) = DDkappa1.block<3, 1>(4, 3).transpose();
    DDkappa1.block<3, 1>(8, 3) = D2kappa1DfDthetae;
    DDkappa1.block<1, 3>(3, 8) = DDkappa1.block<3, 1>(8, 3).transpose();
    DDkappa1.block<3, 1>(0, 7) = -D2kappa1DeDthetaf;
    DDkappa1.block<1, 3>(7, 0) = DDkappa1.block<3, 1>(0, 7).transpose();
    DDkappa1.block<3, 1>(4, 7) = D2kappa1DeDthetaf - D2kappa1DfDthetaf;
    DDkappa1.block<1, 3>(7, 4) = DDkappa1.block<3, 1>(4, 7).transpose();
    DDkappa1.block<3, 1>(8, 7) = D2kappa1DfDthetaf;
    DDkappa1.block<1, 3>(7, 8) = DDkappa1.block<3, 1>(8, 7).transpose();

    assert(isSymmetric(DDkappa1));

    DDkappa2.block<3, 3>(0, 0) = D2kappa2De2;
    DDkappa2.block<3, 3>(0, 4) = -D2kappa2De2 + D2kappa2DeDf;
    DDkappa2.block<3, 3>(4, 0) = -D2kappa2De2 + D2kappa2DfDe;
    DDkappa2.block<3, 3>(4, 4) =
        D2kappa2De2 - (D2kappa2DeDf + D2kappa2DfDe) + D2kappa2Df2;
    DDkappa2.block<3, 3>(0, 8) = -D2kappa2DeDf;
    DDkappa2.block<3, 3>(8, 0) = -D2kappa2DfDe;
    DDkappa2.block<3, 3>(4, 8) = D2kappa2DeDf - D2kappa2Df2;
    DDkappa2.block<3, 3>(8, 4) = D2kappa2DfDe - D2kappa2Df2;
    DDkappa2.block<3, 3>(8, 8) = D2kappa2Df2;
    DDkappa2(3, 3) = D2kappa2Dthetae2;
    DDkappa2(7, 7) = D2kappa2Dthetaf2;
    DDkappa2(3, 7) = DDkappa2(7, 3) = 0.;
    DDkappa2.block<3, 1>(0, 3) = -D2kappa2DeDthetae;
    DDkappa2.block<1, 3>(3, 0) = DDkappa2.block<3, 1>(0, 3).transpose();
    DDkappa2.block<3, 1>(4, 3) = D2kappa2DeDthetae - D2kappa2DfDthetae;
    DDkappa2.block<1, 3>(3, 4) = DDkappa2.block<3, 1>(4, 3).transpose();
    DDkappa2.block<3, 1>(8, 3) = D2kappa2DfDthetae;
    DDkappa2.block<1, 3>(3, 8) = DDkappa2.block<3, 1>(8, 3).transpose();
    DDkappa2.block<3, 1>(0, 7) = -D2kappa2DeDthetaf;
    DDkappa2.block<1, 3>(7, 0) = DDkappa2.block<3, 1>(0, 7).transpose();
    DDkappa2.block<3, 1>(4, 7) = D2kappa2DeDthetaf - D2kappa2DfDthetaf;
    DDkappa2.block<1, 3>(7, 4) = DDkappa2.block<3, 1>(4, 7).transpose();
    DDkappa2.block<3, 1>(8, 7) = D2kappa2DfDthetaf;
    DDkappa2.block<1, 3>(7, 8) = DDkappa2.block<3, 1>(8, 7).transpose();

    assert(isSymmetric(DDkappa2));

    DDkappa3.block<3, 3>(0, 0) = D2kappa3De2;
    DDkappa3.block<3, 3>(0, 4) = -D2kappa3De2 + D2kappa3DeDf;
    DDkappa3.block<3, 3>(4, 0) = -D2kappa3De2 + D2kappa3DfDe;
    DDkappa3.block<3, 3>(4, 4) =
        D2kappa3De2 - (D2kappa3DeDf + D2kappa3DfDe) + D2kappa3Df2;
    DDkappa3.block<3, 3>(0, 8) = -D2kappa3DeDf;
    DDkappa3.block<3, 3>(8, 0) = -D2kappa3DfDe;
    DDkappa3.block<3, 3>(4, 8) = D2kappa3DeDf - D2kappa3Df2;
    DDkappa3.block<3, 3>(8, 4) = D2kappa3DfDe - D2kappa3Df2;
    DDkappa3.block<3, 3>(8, 8) = D2kappa3Df2;
    DDkappa3(3, 3) = D2kappa3Dthetae2;
    DDkappa3(7, 7) = D2kappa3Dthetaf2;
    DDkappa3(3, 7) = DDkappa3(7, 3) = 0.;
    DDkappa3.block<3, 1>(0, 3) = -D2kappa3DeDthetae;
    DDkappa3.block<1, 3>(3, 0) = DDkappa3.block<3, 1>(0, 3).transpose();
    DDkappa3.block<3, 1>(4, 3) = D2kappa3DeDthetae - D2kappa3DfDthetae;
    DDkappa3.block<1, 3>(3, 4) = DDkappa3.block<3, 1>(4, 3).transpose();
    DDkappa3.block<3, 1>(8, 3) = D2kappa3DfDthetae;
    DDkappa3.block<1, 3>(3, 8) = DDkappa3.block<3, 1>(8, 3).transpose();
    DDkappa3.block<3, 1>(0, 7) = -D2kappa3DeDthetaf;
    DDkappa3.block<1, 3>(7, 0) = DDkappa3.block<3, 1>(0, 7).transpose();
    DDkappa3.block<3, 1>(4, 7) = D2kappa3DeDthetaf - D2kappa3DfDthetaf;
    DDkappa3.block<1, 3>(7, 4) = DDkappa3.block<3, 1>(4, 7).transpose();
    DDkappa3.block<3, 1>(8, 7) = D2kappa3DfDthetaf;
    DDkappa3.block<1, 3>(7, 8) = DDkappa3.block<3, 1>(8, 7).transpose();

    assert(isSymmetric(DDkappa3));
  }

  setDependentsDirty();
}
}  // namespace strandsim
