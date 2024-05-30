/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "VertexFaceCollision.hh"

#include "../Core/ElasticStrand.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Utils/Distances.hh"
#include "../Utils/MathUtilities.hh"
#include "../Utils/TextLog.hh"
#include "CollisionUtils.hh"
#include "ElementProxy.hh"

namespace strandsim {

static const double SQ_TOLERANCE = 1e-12;
static const double EXTRA_RADIUS = 1e-3;

void VertexFaceCollision::print(std::ostream& os) const {
  os << "VertexFaceCollision: strand vertex " << m_firstStrand->getGlobalIndex()
     << ' ' << m_firstVertex << " vs. " << m_faceProxy << " by "
     << -m_normalRelativeDisplacement << " at time = " << m_time;

  os << '\n';
  os << "Normal displacement = " << m_normal << '\n';
  os << "Vertex moved from: "
     << m_firstStrand->getFutureState().getVertex(m_firstVertex) << " to "
     << m_firstStrand->getVertex(m_firstVertex);
  os << '\n';
  os << "Face moved: " << *m_faceProxy << '\n';
  os << "Normal relative displacement: " << m_normalRelativeDisplacement;
  os << '\n';
}

bool VertexFaceCollision::analyse() {
  double times[4];

  const int num_verts = m_firstStrand->getNumVertices();
  const Vec3x p_off = m_firstStrand->getVertex(m_firstVertex);
  // NB after taking off the offset p = 0
  Vec3x pf0 = m_faceProxy->getVertex(0) - p_off;
  Vec3x pf1 = m_faceProxy->getVertex(1) - p_off;
  Vec3x pf2 = m_faceProxy->getVertex(2) - p_off;

  const Vec3x dp = m_firstStrand->dynamics().getDisplacement(m_firstVertex);
  Vec3x df0 = m_faceProxy->getDisplacement(0);
  Vec3x df1 = m_faceProxy->getDisplacement(1);
  Vec3x df2 = m_faceProxy->getDisplacement(2);

  const Vec3x dfc = (df0 + df1 + df2) / 3.;

  int fnsign =
      m_faceProxy->knowsNormalSign(true, *m_firstStrand, m_firstVertex);
  m_normal = m_faceProxy->getNormal();

  const Scalar colRadius =
      m_firstStrand->collisionParameters().externalCollisionsRadius(
          m_firstVertex, M_PI / 2);
  const Vec3x dpc = dp - dfc;
  const Scalar disp_rel = dpc.dot(m_normal);
  const Scalar tang_rel =
      sqrt(std::max(0., dpc.squaredNorm() - disp_rel * disp_rel));

  Scalar extraRadius;

  const Scalar contactAngle = m_firstStrand->getContactAngle();
  extraRadius =
      colRadius + (1. + 0.5 * contactAngle) *
                      sqrt(m_firstStrand->getCurrentFlowDOFArea(m_firstVertex));

  m_adhesive_force = 0.;
  m_yield = 0.;
  m_eta = 0.;
  m_power = 1.;

  if (extraRadius < 0.) return false;

  //    const Vec3x& prox = extraRadius * fnsign * m_normal;
  //    pf0 += prox;
  //    pf1 += prox;
  //    pf2 += prox;
  //    df0 += prox;
  //    df1 += prox;
  //    df2 += prox;

  unsigned num_times;
  getCoplanarityTimes(-dp, pf0 - df0, pf1 - df1, pf2 - df2, Vec3x(), pf0, pf1,
                      pf2, times, NULL, num_times);

  for (size_t j = 0; j < num_times; ++j) {
    // TODO: Use barycentric coordinates or point-triangle closest point <
    // epsilon here? closest point < epsilon really just extends the triangle a
    // little bit. Determine if the collision actually happens
    const Scalar dtime = times[j] - 1.0;
    const Vec3x pcol = dtime * (dp);
    const Vec3x f0col = pf0 + dtime * df0;
    const Vec3x f1col = pf1 + dtime * df1;
    const Vec3x f2col = pf2 + dtime * df2;

    const Vec3x cp = ClosestPtPointTriangle(pcol, f0col, f1col, f2col);

    const Scalar pqdist2 = (pcol - cp).squaredNorm();
    // If, when they are coplanar, the objects are sufficiently close, register
    // a collision
    if (pqdist2 < extraRadius * extraRadius) {
      m_time = times[j];

      const Scalar cpA = m_firstStrand->getCurrentFlowDOFArea(m_firstVertex);

      if (cpA < 1e-12 || disp_rel < -1e-7) {
        m_do_soc_solve = pqdist2 < colRadius * colRadius;
      } else {
        m_do_soc_solve = true;
      }

      //            Scalar m_u, m_v, m_w;
      computeBarycentricCoordinates(f0col, f1col, f2col, pcol, m_u, m_v, m_w);
      // computeBarycentricCoordinates coords could be outside of [0,1] right
      // now because we've extended the triangles a little bit
      assert(isSmall(m_u + m_v + m_w - 1.0));

      m_meshDisplacement = m_u * df0 + m_v * df1 + m_w * df2;
      const Vec3x relativeDisplacement =
          (1 - m_time) * (dp - m_meshDisplacement);

      //            m_offset = relativeDisplacement ;
      m_offset = (m_u * pf0 + m_v * pf1 + m_w * pf2) -
                 m_meshDisplacement  // orig point on mesh
                 + dp;               // orig point on rod

      const Scalar nDepl = relativeDisplacement.dot(m_normal);
      if (!fnsign) {
        // Normal sign was unknown, now we know that it should be opposed to
        // relativeDisplacement
        fnsign = (nDepl > 0. ? -1 : 1);
        m_faceProxy->setNormalSign(fnsign, m_time, *m_firstStrand,
                                   m_firstVertex);
        //                m_meshDisplacement += extraRadius * fnsign * m_normal
        //                ;
      } else {
        if (fnsign * nDepl > 0.) {
          return false;
        }

        //                const Scalar belowRadius = m_offset.dot( fnsign *
        //                m_normal ) + extraRadius ; if( belowRadius > 0. &&
        //                extraRadius > 0. )
        //                {
        //                    // We're already closer than the collision radius,
        //                    let proximity handle that m_meshDisplacement -=
        //                    prox ; m_offset.setZero() ;
        //                }
      }
      m_normal = fnsign * m_normal;

      Scalar l0 = 0.;
      Scalar l1 = 0.;

      m_adhesive_force = 0.;
      m_power = 1.;
      m_yield = 0.;
      m_eta = 0.;

      const Scalar theta = m_firstStrand->getContactAngle();
      const Scalar cpr = sqrt(m_firstStrand->getRadiusA(m_firstVertex) *
                              m_firstStrand->getRadiusB(m_firstVertex));

      if (cpA > 1e-12) {
        // test prev edge center
        VecXx color = m_firstStrand->getFlowNewComponents(m_firstVertex);

        const Scalar eta = m_firstStrand->getFlowConsistencyIndex(color);
        const Scalar n = m_firstStrand->getFlowBehaviorIndex(color);
        const Scalar tilde_sigma_Y =
            m_firstStrand->getFlowYieldStress(color) * 0.8164965809;

        const Scalar relative_vel = disp_rel / m_firstStrand->getDt();
        const Scalar tang_vel = tang_rel / m_firstStrand->getDt();

        m_power = n;

        if (m_firstVertex > 0) {
          int m_prevVertex = m_firstVertex - 1;
          const Vec3x p_prev = m_firstStrand->getVertex(m_prevVertex) - p_off;
          const Vec3x dp_prev =
              m_firstStrand->dynamics().getDisplacement(m_prevVertex);
          Vec3x pcol_prev = (p_prev + dtime * dp_prev + pcol) * 0.5;
          Scalar A0 =
              (m_firstStrand->getCurrentFlowDOFArea(m_prevVertex) + cpA) * 0.5;

          const Vec3x cp_prev =
              ClosestPtPointTriangle(pcol_prev, f0col, f1col, f2col);
          const Scalar d0 =
              std::max(0., (pcol_prev - cp_prev).dot(m_normal)) + extraRadius;
          const Scalar l = (cp_prev - cp).norm();

          const Scalar abs0 =
              wedge_wet_abscissa(cpA, A0, extraRadius, d0, theta);

          const Scalar Ac = cpA * (1. - abs0 * 0.5) + A0 * abs0 * 0.5;
          const Scalar dc = extraRadius * (1. - abs0 * 0.5) + d0 * abs0 * 0.5;

          const Scalar diff_vel =
              sgn(relative_vel) *
              pow(fabs(relative_vel) / std::max(cpr, dc), n);

          const Scalar elastic_force =
              0.75 * M_PI * cpr * (eta * diff_vel + tilde_sigma_Y);

          l0 = abs0 * l;

          m_adhesive_force +=
              (m_firstStrand->collisionParameters().adhesionForcePlanar(Ac, dc,
                                                                        color) +
               elastic_force) *
              l0;
          m_yield += 0.75 * M_PI * cpr * l0 * tilde_sigma_Y;
          m_eta += 0.75 * M_PI * cpr * eta / pow(std::max(cpr, dc), n) * l0;
        }

        // text next edge center
        if (m_firstVertex < num_verts - 1) {
          int m_nextVertex = m_firstVertex + 1;
          const Vec3x p_next = m_firstStrand->getVertex(m_nextVertex) - p_off;
          const Vec3x dp_next =
              m_firstStrand->dynamics().getDisplacement(m_nextVertex);
          Vec3x pcol_next = (p_next + dtime * dp_next + pcol) * 0.5;
          Scalar A1 =
              (m_firstStrand->getCurrentFlowDOFArea(m_nextVertex) + cpA) * 0.5;

          const Vec3x cp_next =
              ClosestPtPointTriangle(pcol_next, f0col, f1col, f2col);
          const Scalar d1 =
              std::max(0., (pcol_next - cp_next).dot(m_normal)) + extraRadius;
          const Scalar l = (cp_next - cp).norm();

          const Scalar abs1 =
              wedge_wet_abscissa(cpA, A1, extraRadius, d1, theta);

          const Scalar Ac = cpA * (1. - abs1 * 0.5) + A1 * abs1 * 0.5;
          const Scalar dc = extraRadius * (1. - abs1 * 0.5) + d1 * abs1 * 0.5;

          const Scalar diff_vel =
              sgn(relative_vel) *
              pow(fabs(relative_vel) / std::max(cpr, dc), n);

          const Scalar elastic_force =
              0.75 * M_PI * cpr * (eta * diff_vel + tilde_sigma_Y);

          l1 = abs1 * l;

          m_adhesive_force +=
              (m_firstStrand->collisionParameters().adhesionForcePlanar(Ac, dc,
                                                                        color) +
               elastic_force) *
              l1;
          m_yield += 0.75 * M_PI * cpr * l1 * tilde_sigma_Y;
          m_eta += 0.75 * M_PI * cpr * eta / pow(std::max(cpr, dc), n) * l1;
        }
      }

      postAnalyse(relativeDisplacement);

      return true;
    }
  }
  return false;
}

bool compare(const VertexFaceCollision* vf1, const VertexFaceCollision* vf2) {
  if (vf1->m_firstStrand == vf2->m_firstStrand)
    if (vf1->m_firstVertex == vf2->m_firstVertex)
      if (vf1->m_faceProxy == vf2->m_faceProxy)
        return vf1->m_time < vf2->m_time;
      else
        return vf1->m_faceProxy < vf2->m_faceProxy;
    else
      return vf1->m_firstVertex < vf2->m_firstVertex;
  else
    return vf1->m_firstStrand < vf2->m_firstStrand;
}

}  // namespace strandsim
