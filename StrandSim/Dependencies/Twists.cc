/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Twists.hh"

#include "../Core/ElasticStrandUtils.hh"
#include "../Utils/MathUtilities.hh"

//#define CURVATURE_USE_SINTHETA

#ifndef CURVATURE_USE_SINTHETA
#define CURVATURE_USE_TANTHETA
#endif

namespace strandsim {

void Twists::compute() {
  m_value.resize(m_size);
  const std::vector<Scalar>& refTwists = m_refTwists.get();
  const VecXx& dofs = m_dofs.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    m_value[vtx] = refTwists[vtx] + dofs[4 * vtx + 3] - dofs[4 * vtx - 1];
  }

  setDependentsDirty();
}

void GradTwists::compute() {
  m_value.resize(m_size);
  const Vec3xArray& curvatureBinormals = m_curvatureBinormals.get();
  const std::vector<Scalar>& lengths = m_lengths.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    Vec11x& Dtwist = m_value[vtx];

    const Vec3x& kb = curvatureBinormals[vtx];

    Dtwist.segment<3>(0) = -0.5 / lengths[vtx - 1] * kb;
    Dtwist.segment<3>(8) = 0.5 / lengths[vtx] * kb;
    Dtwist.segment<3>(4) = -(Dtwist.segment<3>(0) + Dtwist.segment<3>(8));
    Dtwist(3) = -1;
    Dtwist(7) = 1;
  }

  setDependentsDirty();
}

void GradTwistsSquared::compute() {
  m_value.resize(m_size);
  const Vec11xArray& gradTwists = m_gradTwists.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    const Vec11x& gradTwist = gradTwists[vtx];
    m_value[vtx] = gradTwist * gradTwist.transpose();
  }

  setDependentsDirty();
}

void HessTwists::compute() {
  m_value.resize(m_size);
  const Vec3xArray& tangents = m_tangents.get();
  const std::vector<Scalar>& lengths = m_lengths.get();
  const Vec3xArray& curvatureBinormals = m_curvatureBinormals.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    Mat11x& DDtwist = m_value[vtx];

    const Vec3x& te = tangents[vtx - 1];
    const Vec3x& tf = tangents[vtx];
    const Scalar norm_e = lengths[vtx - 1];
    const Scalar norm_f = lengths[vtx];
    const Vec3x& kb = curvatureBinormals[vtx];
#ifdef CURVATURE_USE_SINTHETA
    Scalar chi = (te + tf).norm();

    //    assert( chi>0 );
    if (chi < 1e-12) {
      WarningStream(g_log, "")
          << "StrandState::computeHessTwist for state " << this
          << " chi = " << chi << " te = " << te << " tf = " << tf;
      chi = 1e-12;
    }

    const Vec3x& tilde_t = 1.0 / chi * (te + tf);

    const Mat3x& D2mDe2 =
        -0.5 / (chi * square(norm_e)) *
        symPart(outerProd<3>(kb, tilde_t + (1.0 + tilde_t.dot(tf)) * te));
    const Mat3x& D2mDf2 =
        -0.5 / (chi * square(norm_f)) *
        symPart(outerProd<3>(kb, tilde_t + (1.0 + tilde_t.dot(te)) * tf));
    const Mat3x& D2mDeDf =
        1.0 / (norm_e * norm_f * chi) *
        (crossMat(te) - 0.5 * outerProd<3>(kb, tilde_t + tilde_t.dot(te) * tf));
    const Mat3x& D2mDfDe = D2mDeDf.transpose();
#else
    Scalar chi = 1 + te.dot(tf);

    //    assert( chi>0 );
    if (chi < 1e-12) {
      WarningStream(g_log, "")
          << "StrandState::computeHessTwist for state " << this
          << " chi = " << chi << " te = " << te << " tf = " << tf;
      chi = 1e-12;
    }

    const Vec3x& tilde_t = 1.0 / chi * (te + tf);

    const Mat3x& D2mDe2 =
        -0.25 / square(norm_e) *
        (outerProd<3>(kb, te + tilde_t) + outerProd<3>(te + tilde_t, kb));
    const Mat3x& D2mDf2 =
        -0.25 / square(norm_f) *
        (outerProd<3>(kb, tf + tilde_t) + outerProd<3>(tf + tilde_t, kb));
    const Mat3x& D2mDeDf =
        0.5 / (norm_e * norm_f) *
        (2.0 / chi * crossMat(te) - outerProd<3>(kb, tilde_t));
    const Mat3x& D2mDfDe = D2mDeDf.transpose();
#endif
    DDtwist.block<3, 3>(0, 0) = D2mDe2;
    DDtwist.block<3, 3>(0, 4) = -D2mDe2 + D2mDeDf;
    DDtwist.block<3, 3>(4, 0) = -D2mDe2 + D2mDfDe;
    DDtwist.block<3, 3>(4, 4) = D2mDe2 - (D2mDeDf + D2mDfDe) + D2mDf2;
    DDtwist.block<3, 3>(0, 8) = -D2mDeDf;
    DDtwist.block<3, 3>(8, 0) = -D2mDfDe;
    DDtwist.block<3, 3>(8, 4) = D2mDfDe - D2mDf2;
    DDtwist.block<3, 3>(4, 8) = D2mDeDf - D2mDf2;
    DDtwist.block<3, 3>(8, 8) = D2mDf2;

    assert(isSymmetric(DDtwist));
  }

  setDependentsDirty();
}

}  // namespace strandsim
