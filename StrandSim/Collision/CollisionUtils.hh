/**
 * \copyright 2010 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSIM_COLLISIONUTILS_HH
#define STRANDSIM_COLLISIONUTILS_HH

#include "../Core/CollisionParameters.hh"
#include "../Core/Definitions.hh"

namespace strandsim {

class ElasticStrand;
class StrandState;
class ProximityCollisionDatabase;

/////////////////////////////////////////////////////////////////
// Collision detection code adapted from Robert Bridson's website

void addUnique(std::vector<double>& a, double e);

double triple(const Vec3x& a, const Vec3x& b, const Vec3x& c);

double signed_volume(const Vec3x& x0, const Vec3x& x1, const Vec3x& x2,
                     const Vec3x& x3);

void getCoplanarityTimes(const Vec3x& x0, const Vec3x& x1, const Vec3x& x2,
                         const Vec3x& x3, const Vec3x& xnew0,
                         const Vec3x& xnew1, const Vec3x& xnew2,
                         const Vec3x& xnew3, double* times, double* errors,
                         unsigned& num_times);

void getIntersectionPoint(const Vec3x& x0, const Vec3x& xnew0,
                          const Vec3x& xnew1, const Vec3x& xnew2,
                          const Vec3x& xnew3, double* times, double* errors,
                          unsigned& num_times);

void buildFrame(const Vec3x& n_hat, Vec3x& t1, Vec3x& t2);

/////
// return min_( t in 0, 2 Pi ) | E(t) dot normal |
// where E(t) describes the bounding ellipsoid of the strand at edge
Scalar getEllipticExternalCollisionOffset(const ElasticStrand& strand,
                                          const unsigned edge,
                                          const Vec3x& normal);

bool analyseRoughRodRodCollision(const ProximityCollisionDatabase& database,
                                 const ElasticStrand* sP,
                                 const ElasticStrand* sQ, const int iP,
                                 const int iQ, const Scalar& contactAngle,
                                 Vec3x& normalQtoP, Scalar& s, Scalar& t,
                                 Scalar& distance, Scalar& adhesion,
                                 Scalar& yield, Scalar& eta, Scalar& power,
                                 Scalar& relative_vel, bool& solidTouching);

}  // namespace strandsim

#endif

// TODO:
//   o It would be nice to handle degenerate cases better in these methods.
//     They all handle degenerate cases, but getting PREDICTABLE behavior out
//     would rock!!!
