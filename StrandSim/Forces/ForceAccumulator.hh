/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef FORCEACCUMULATOR_HH_
#define FORCEACCUMULATOR_HH_

#include "../Core/BandMatrix.hh"
#include "../Core/Definitions.hh"

namespace strandsim {

class ElasticStrand;
class StrandState;

template <typename ForceT>
class ForceAccumulator {
 public:
  template <typename AccumulatedT>
  static void accumulateCurrent(AccumulatedT& accumulated,
                                ElasticStrand& strand) {
    // Yes, this statement does nothing... It is there to guarantee that
    // accumulateCurrent will never be called on dissipative forces, as they
    // require evolution to be estimated. So if you get a compilation error
    // here, that's why. If thas is a problem just comment out the line below,
    // but make sure that dissipative forces still make sense when called on the
    // current state (most likely they should return zero).
    ForceT::NonDissipativeForce;

    accumulate(accumulated, strand, strand.getCurrentState());
  }

  static void accumulateMultipliers(VecXx& multipliers,
                                    const ElasticStrand& strand) {
    typename ForceT::LocalMultiplierType localL;
    for (IndexType vtx = ForceT::s_first;
         vtx < strand.getNumVertices() - ForceT::s_last; ++vtx) {
      ForceT::computeLocal(localL, strand, vtx);
      ForceT::addInPositionMultiplier(multipliers, vtx, localL);
    }
  }

  template <typename AccumulatedT>
  static void accumulateFuture(AccumulatedT& accumulated,
                               ElasticStrand& strand) {
    accumulate(accumulated, strand, strand.getFutureState());
  }

  static void accumulate(Scalar& energy, const ElasticStrand& strand,
                         StrandState& state) {
    //        assert( state.m_framesUpToDate );

    for (IndexType vtx = ForceT::s_first;
         vtx < state.m_numVertices - ForceT::s_last; ++vtx) {
      energy += ForceT::localEnergy(strand, state, vtx);
    }
  }

  static void accumulate(VecXx& force, const ElasticStrand& strand,
                         StrandState& state) {
    //        assert( state.m_framesUpToDate );

    typename ForceT::LocalForceType localF;
    for (IndexType vtx = ForceT::s_first;
         vtx < state.m_numVertices - ForceT::s_last; ++vtx) {
      ForceT::computeLocal(localF, strand, state, vtx);
      ForceT::addInPosition(force, vtx, localF);
    }
  }

  static void accumulate(JacobianMatrixType& Jacobian,
                         const ElasticStrand& strand, StrandState& state) {
    //        assert( state.m_framesUpToDate );

    typename ForceT::LocalJacobianType localJ;
    for (IndexType vtx = ForceT::s_first;
         vtx < state.m_numVertices - ForceT::s_last; ++vtx) {
      ForceT::computeLocal(localJ, strand, state, vtx);
      ForceT::addInPosition(Jacobian, vtx, localJ);
    }
  }
};

}  // namespace strandsim

#endif /* FORCEACCUMULATOR_HH_ */
