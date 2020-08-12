/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef VISCOUSORNOTVISCOUS_HH_
#define VISCOUSORNOTVISCOUS_HH_

#include "../Core/ElasticStrand.hh"

namespace strandsim {

// These classes are taken as template arguments for the internal forces,
// indicating whether we want the non-viscous or the viscous version.
// The forces call their ViscousT's static methods returning the appropriate
// stiffness and "rest shape" (the actual rest-shape for non-viscous or the
// shape at the beginning of time step for viscous).

class NonViscous {
 protected:
  NonViscous() {}

  virtual ~NonViscous() {}

 public:
  static std::string getName() { return ""; }

  static Scalar bendingCoefficient(const ElasticStrand& strand, int vtx) {
    return strand.m_parameters.bendingCoefficient(vtx, strand.getNumVertices());
  }
  static Mat2x bendingMatrix(const ElasticStrand& strand, int vtx) {
    return strand.m_parameters.bendingMatrix(vtx, strand.getNumVertices());
  }

  static const Vec4x kappaBar(const ElasticStrand& strand, int vtx) {
    return strand.m_restKappas[vtx];
  }

  static Scalar kt(const ElasticStrand& strand, int vtx) {
    return strand.m_parameters.getKt(vtx, strand.getNumVertices());
  }

  static Scalar thetaBar(const ElasticStrand& strand, int vtx) {
    return strand.m_restTwists[vtx];
  }

  static Scalar ks(const ElasticStrand& strand, int vtx) {
    return strand.m_parameters.getKs(vtx, strand.getNumVertices());
  }

  static Scalar kmb(const ElasticStrand& strand, int vtx) {
    return strand.m_parameters.getKs(vtx, strand.getNumVertices()) *
           strand.m_parameters.getMinBendingMultiplier();
  }

  static Scalar kf(const ElasticStrand& strand, int vtx) {
    return strand.m_parameters.getKf(vtx, strand.getNumVertices());
  }

  static Scalar ellBar(const ElasticStrand& strand, int vtx) {
    return strand.m_restLengths[vtx];
  }

  static Scalar neighborDist(const ElasticStrand& strand, int vtx) {
    return strand.m_restNeighborDistances[vtx];
  }

  static bool isViscous() { return false; }

  class NonDissipativeForce {};
};

class Viscous {
 protected:
  Viscous() {}

  virtual ~Viscous() {}

 public:
  static std::string getName() { return "viscous "; }

  static Scalar bendingCoefficient(const ElasticStrand& strand, int vtx) {
    return strand.m_parameters.viscousBendingCoefficient(
        vtx, strand.getNumVertices());
  }

  static Mat2x bendingMatrix(const ElasticStrand& strand, int vtx) {
    return strand.m_parameters.viscousBendingMatrix(vtx,
                                                    strand.getNumVertices());
  }

  static const Vec4x kappaBar(const ElasticStrand& strand, int vtx) {
    return strand.m_currentState->m_kappas[vtx];
  }

  static Scalar kt(const ElasticStrand& strand, int vtx) {
    return strand.m_parameters.getViscousKt(vtx, strand.getNumVertices());
  }

  static Scalar thetaBar(const ElasticStrand& strand, int vtx) {
    return strand.m_currentState->m_twists[vtx];
  }

  static Scalar ks(const ElasticStrand& strand, int vtx) {
    return strand.m_parameters.getViscousKs(vtx, strand.getNumVertices());
  }

  static Scalar kmb(const ElasticStrand& strand, int vtx) {
    return strand.m_parameters.getViscousKs(vtx, strand.getNumVertices()) *
           strand.m_parameters.getMinBendingMultiplier();
  }

  static Scalar kf(const ElasticStrand& strand, int vtx) {
    return strand.m_parameters.getViscousKf(vtx, strand.getNumVertices());
  }

  static Scalar ellBar(const ElasticStrand& strand, int vtx) {
    return strand.m_currentState->m_lengths[vtx];
  }

  static Scalar neighborDist(const ElasticStrand& strand, int vtx) {
    return strand.m_currentState->m_neighbor_distances[vtx];
  }

  static bool isViscous() { return true; }

  class DissipativeForce {};
};

}  // namespace strandsim

#endif /* VISCOUSORNOTVISCOUS_HH_ */
