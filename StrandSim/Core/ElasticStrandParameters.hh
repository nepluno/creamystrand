/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ELASTICSTRANDPARAMETERS_HH_
#define ELASTICSTRANDPARAMETERS_HH_

#include <boost/serialization/access.hpp>
#include <boost/serialization/version.hpp>

#include "../Dependencies/BendingProducts.hh"
#include "Definitions.hh"

namespace strandsim {

class ElasticStrandParameters {
 public:
  ElasticStrandParameters(int paramIndex = -1);

  ElasticStrandParameters(const ElasticStrandParameters& other);

  ElasticStrandParameters(const ElasticStrandParameters& other,
                          const VecXx& radiusMultiplier);

  ElasticStrandParameters(int paramIndex, Scalar radiusA, Scalar radiusB,
                          Scalar YoungsModulus, Scalar shearModulus,
                          Scalar density, Scalar viscosity, Scalar airDrag,
                          Scalar baseRotation, Scalar stretchMultiplier,
                          Scalar fixingMultiplier, Scalar minBendingAngle,
                          Scalar maxFlowGradientRatio,
                          Scalar minBendingMultiplier,
                          bool flowOpenAtStart = true,
                          bool flowOpenAtEnd = true);

  ElasticStrandParameters(
      int paramIndex, Scalar radiusA, Scalar radiusB, Scalar YoungsModulus,
      Scalar shearModulus, Scalar density, Scalar viscosity, Scalar airDrag,
      Scalar baseRotation, Scalar stretchMultiplier, Scalar fixingMultiplier,
      Scalar minBendingAngle, Scalar maxFlowGradientRatio,
      const VecXx& radiusMultiplier, Scalar minBendingMultiplier,
      bool flowOpenAtStart = true, bool flowOpenAtEnd = true);

  ElasticStrandParameters& operator=(const ElasticStrandParameters& other);

  Scalar getBaseRotation() const { return m_baseRotation.get(); }
  Scalar getDensity() const { return m_density; }
  Scalar getMaxRadiusMultiplier() const { return m_maxRadiusMultiplier; }
  Scalar getKf(int vtx, int num_verts) const {
    return getKs(vtx, num_verts) * m_fixingMultiplier;
  }
  Scalar getKs(int vtx, int num_verts) const {
    const Scalar interpol = m_internalRadiusMultiplier *
                            interpolatedRadiusMultiplier(vtx, num_verts);

    return interpol * interpol * m_ks.get() * m_stretchMultiplier;
  }
  Scalar getKt(int vtx, int num_verts) const {
    const Scalar interpol = m_internalRadiusMultiplier *
                            interpolatedRadiusMultiplier(vtx, num_verts);
    const Scalar interpol2 = interpol * interpol;

    return interpol2 * interpol2 * m_kt.get();
  }
  Scalar getRadiusA(int vtx, int num_verts) const {
    return m_internalRadiusMultiplier *
           interpolatedRadiusMultiplier(vtx, num_verts) *
           m_ellipticalRadii.get().first;
  }
  Scalar getRadiusB(int vtx, int num_verts) const {
    return m_internalRadiusMultiplier *
           interpolatedRadiusMultiplier(vtx, num_verts) *
           m_ellipticalRadii.get().second;
  }
  Scalar getRadiusAtAngle(int vtx, int num_verts, Scalar angle) const {
    return getRadiusB(vtx, num_verts) * getRadiusShrinkAtAngle(angle);
  }
  Scalar getShearModulus() const { return m_shearModulus.get(); }
  Scalar getYoungsModulus() const { return m_youngsModulus.get(); }
  Scalar getViscosity() const { return m_viscosity; }
  Scalar getAirDrag() const { return m_airDrag; }
  Scalar getMinBendingAngle() const { return m_minBendingAngle; }
  Scalar getMaxFlowGradientRatio() const { return m_maxFlowGradientRatio; }

  Scalar bendingCoefficient(int vtx, int num_verts) const {
    const Scalar interpol = m_internalRadiusMultiplier *
                            interpolatedRadiusMultiplier(vtx, num_verts);
    const Scalar interpol2 = interpol * interpol;

    return interpol2 * interpol2 * m_youngsModulus.get();
  }

  Mat2x bendingMatrix(int vtx, int num_verts) const {
    return bendingCoefficient(vtx, num_verts) * m_bendingMatrixBase.get();
  }

  void setBaseRotation(Scalar baseRotation);
  void setDensity(Scalar density);
  void setRadii(Scalar radiusA, Scalar radiusB);
  void setShearModulus(Scalar shearModulus);
  void setYoungsModulus(Scalar YoungsModulus);
  void setViscosity(Scalar viscosity);
  void setAirDrag(Scalar airDrag);

  //! Since we dont' want to save dirty values, we need a way to compute them
  //! before
  void computeDependencies();

  void computeViscousForceCoefficients(Scalar dt);

  void setRadiusMultiplier(int vtx, double radiusMultiplier) {
    m_radiusMultiplier(vtx) = radiusMultiplier;
  }

  void resizeRadiusMultiplier(int num) {
    m_radiusMultiplier.conservativeResize(num);
  }

  Scalar viscousBendingCoefficient(int vtx, int num_verts) const {
    const Scalar interpol = m_internalRadiusMultiplier *
                            interpolatedRadiusMultiplier(vtx, num_verts);
    const Scalar interpol2 = interpol * interpol;

    return interpol2 * interpol2 * m_viscousBendingCoefficientBase;
  }

  Mat2x viscousBendingMatrix(int vtx, int num_verts) const {
    return viscousBendingCoefficient(vtx, num_verts) *
           m_bendingMatrixBase.get();
  }

  Scalar getViscousKs(int vtx, int num_verts) const {
    const Scalar interpol = m_internalRadiusMultiplier *
                            interpolatedRadiusMultiplier(vtx, num_verts);

    return interpol * interpol * m_viscousKs * m_stretchMultiplier;
  }

  Scalar getViscousKf(int vtx, int num_verts) const {
    return getViscousKs(vtx, num_verts) * m_fixingMultiplier;
  }

  Scalar getViscousKt(int vtx, int num_verts) const {
    const Scalar interpol = m_internalRadiusMultiplier *
                            interpolatedRadiusMultiplier(vtx, num_verts);
    const Scalar interpol2 = interpol * interpol;

    return interpol2 * interpol2 * m_viscousKt;
  }

  Scalar getSlipLength() const { return m_slipLength; }

  Scalar getViscousLamb() const { return m_viscousLamb; }

  Scalar getDt() const { return m_dt; }

  Scalar getKs() const { return m_ks.get(); }

  Scalar getKf() const { return getKs() * m_fixingMultiplier; }

  Scalar getRawRadius() const { return m_ellipticalRadii.get().first; }

  void setKs(const Scalar ks) { m_ks.set(ks); }

  BendingMatrixBase& getBendingMatrixBase() { return m_bendingMatrixBase; }

  const Mat2x& bendingMatrixBase() const { return m_bendingMatrixBase.get(); }

  Scalar& dt() { return m_dt; }

  const Scalar& internalRadiusMultiplier() const {
    return m_internalRadiusMultiplier;
  }
  Scalar& internalRadiusMultiplier() { return m_internalRadiusMultiplier; }
  Scalar interpolatedRadiusMultiplier(int vtx, int num_verts) const;

  bool collisionFree() const { return m_collisionFree; }

  void setCollisionFree(bool cf) { m_collisionFree = cf; }

  bool isFlowOpenAtStart() const { return m_flowOpenAtStart; }
  bool isFlowOpenAtEnd() const { return m_flowOpenAtEnd; }

  int getParamIndex() const { return m_paramIndex; }

  Scalar getFixingMultiplier() const { return m_fixingMultiplier; }

  Scalar getMinBendingMultiplier() const { return m_minBendingMultiplier; }

 private:
  friend class CollisionParameters;

  Scalar getRadiusShrinkAtAngle(const Scalar angle) const;

  static Scalar radiusShrinkAtAngle(const Scalar ratio, const Scalar angle);

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& m_ellipticalRadii;
    ar& m_youngsModulus;
    ar& m_shearModulus;
    ar& m_density;
    ar& m_baseRotation;

    // Versioning -- increase BOOST_CLASS_VERSION below each time new properties
    // are added
    if (version > 0) {
      ar& m_viscosity;
      if (version > 1) {
        ar& m_airDrag;
      }
      if (version > 2) {
        ar& m_radiusMultiplier(0);
        ar& m_radiusMultiplier(m_radiusMultiplier.size() - 1);
        ar& m_bendingMatrixBase;
      }
      if (version > 3) {
        ar& m_internalRadiusMultiplier;
      }
    }

    ar& m_ks;
    ar& m_kt;
  }

  double m_density;
  double m_viscosity;
  double m_airDrag;
  double m_slipLength;

  VecXx m_radiusMultiplier;
  double m_maxRadiusMultiplier;

  // Radius multiplier that should not be changed by user input.
  // Used for stiffening rods in reverse hair-do
  // Not copied in operator= or copy constructor
  Scalar m_internalRadiusMultiplier;

  // Used for directly modifying stretch Ks
  Scalar m_stretchMultiplier;

  // Used for directly modifying fixing Ks
  Scalar m_fixingMultiplier;

  Scalar m_minBendingMultiplier;

  // Used for limiting the bending angle
  Scalar m_minBendingAngle;

  // Used for limiting the gradient of flow area (so that the mass ratio between
  // neighbor vertices won't be too large to blow up the simulation)
  Scalar m_maxFlowGradientRatio;

  // Computed viscous force coefficients. Since we are not maintaining m_dt
  // here, we cannot update them automatically whenever a physical parameter is
  // changed. So computeViscousForceCoefficients(dt) MUST be called at the
  // beginning of each time step, in case something has changed.
  double m_viscousBendingCoefficientBase;
  Scalar m_viscousKt;
  Scalar m_viscousKs;
  Scalar m_viscousLamb;  // This makes sense

  // Dependencies
  mutable EllipticalRadii m_ellipticalRadii;
  mutable BaseRotation m_baseRotation;
  mutable BendingMatrixBase m_bendingMatrixBase;
  mutable YoungsModulus m_youngsModulus;
  mutable ShearModulus m_shearModulus;
  mutable ElasticKs m_ks;
  mutable ElasticKt m_kt;

  bool m_collisionFree;

  bool m_flowOpenAtStart;
  bool m_flowOpenAtEnd;

  int m_paramIndex;

  Scalar m_dt;  // This really doesn't belong here
};

}  // namespace strandsim

BOOST_CLASS_VERSION(strandsim::ElasticStrandParameters, 4)

#endif /* ELASTICSTRANDPARAMETERS_HH_ */
