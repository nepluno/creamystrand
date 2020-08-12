/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "ElasticStrandParameters.hh"

#include "ElasticStrandUtils.hh"

namespace strandsim {

ElasticStrandParameters::ElasticStrandParameters(int paramIndex)
    :  // Dummy default values ElasticStrand's method do not generate FPEs
       // when the parameters are not set yet
      m_paramIndex(paramIndex),
      m_density(0.),                   //
      m_viscosity(0.),                 //
      m_airDrag(0.),                   //
      m_internalRadiusMultiplier(1.),  //
      m_ellipticalRadii(0., 0.),       //
      m_baseRotation(0.),              //
      m_minBendingAngle(30.0 * M_PI / 180.0),
      m_bendingMatrixBase(m_ellipticalRadii, m_baseRotation),  //
      m_youngsModulus(0.),                                     //
      m_shearModulus(0.),                                      //
      m_slipLength(0.0002),
      m_stretchMultiplier(1.0),
      m_fixingMultiplier(1.0),
      m_minBendingMultiplier(1.0),
      m_maxFlowGradientRatio(10.0),
      m_ks(m_ellipticalRadii, m_youngsModulus),  //
      m_kt(m_ellipticalRadii, m_shearModulus),
      m_flowOpenAtStart(true),
      m_flowOpenAtEnd(true),
      m_collisionFree(false) {
  m_radiusMultiplier.resize(2);
  m_radiusMultiplier(0) = m_radiusMultiplier(1) = 1.;
  m_maxRadiusMultiplier = m_radiusMultiplier.maxCoeff();
}

ElasticStrandParameters::ElasticStrandParameters(
    const ElasticStrandParameters& other)
    : m_paramIndex(other.m_paramIndex),
      m_density(other.m_density),                    //
      m_viscosity(other.m_viscosity),                //
      m_airDrag(other.m_airDrag),                    //
      m_radiusMultiplier(other.m_radiusMultiplier),  //
      m_internalRadiusMultiplier(1.),                //
      m_ellipticalRadii(other.m_ellipticalRadii.get().first,
                        other.m_ellipticalRadii.get().second),  //
      m_baseRotation(other.m_baseRotation.get()),               //
      m_minBendingAngle(other.m_minBendingAngle),
      m_bendingMatrixBase(m_ellipticalRadii, m_baseRotation),  //
      m_youngsModulus(other.m_youngsModulus.get()),            //
      m_shearModulus(other.m_shearModulus.get()),              //
      m_slipLength(other.m_slipLength),
      m_stretchMultiplier(other.m_stretchMultiplier),
      m_fixingMultiplier(other.m_fixingMultiplier),
      m_minBendingMultiplier(other.m_minBendingMultiplier),
      m_maxFlowGradientRatio(other.m_maxFlowGradientRatio),
      m_ks(m_ellipticalRadii, m_youngsModulus),  //
      m_kt(m_ellipticalRadii, m_shearModulus),
      m_flowOpenAtStart(other.m_flowOpenAtStart),
      m_flowOpenAtEnd(other.m_flowOpenAtEnd),
      m_collisionFree(other.m_collisionFree) {
  m_maxRadiusMultiplier = m_radiusMultiplier.maxCoeff();
}

ElasticStrandParameters::ElasticStrandParameters(
    const ElasticStrandParameters& other, const VecXx& radiusMultiplier)
    : m_paramIndex(other.m_paramIndex),
      m_density(other.m_density),            //
      m_viscosity(other.m_viscosity),        //
      m_airDrag(other.m_airDrag),            //
      m_radiusMultiplier(radiusMultiplier),  //
      m_internalRadiusMultiplier(1.),        //
      m_ellipticalRadii(other.m_ellipticalRadii.get().first,
                        other.m_ellipticalRadii.get().second),  //
      m_baseRotation(other.m_baseRotation.get()),               //
      m_minBendingAngle(other.m_minBendingAngle),
      m_minBendingMultiplier(other.m_minBendingMultiplier),
      m_bendingMatrixBase(m_ellipticalRadii, m_baseRotation),  //
      m_youngsModulus(other.m_youngsModulus.get()),            //
      m_shearModulus(other.m_shearModulus.get()),              //
      m_slipLength(other.m_slipLength),
      m_stretchMultiplier(other.m_stretchMultiplier),
      m_fixingMultiplier(other.m_fixingMultiplier),
      m_maxFlowGradientRatio(other.m_maxFlowGradientRatio),
      m_ks(m_ellipticalRadii, m_youngsModulus),  //
      m_kt(m_ellipticalRadii, m_shearModulus),
      m_flowOpenAtStart(other.m_flowOpenAtStart),
      m_flowOpenAtEnd(other.m_flowOpenAtEnd),
      m_collisionFree(other.m_collisionFree) {
  m_maxRadiusMultiplier = m_radiusMultiplier.maxCoeff();
}

ElasticStrandParameters::ElasticStrandParameters(
    int paramIndex, Scalar radiusA, Scalar radiusB, Scalar YoungsModulus,
    Scalar shearModulus, Scalar density, Scalar viscosity, Scalar airDrag,
    Scalar baseRotation, Scalar stretchMultiplier, Scalar fixingMultiplier,
    Scalar minBendingAngle, Scalar maxFlowGradientRatio,
    Scalar minBendingMultiplier, bool flowOpenAtStart, bool flowOpenAtEnd)
    : m_paramIndex(paramIndex),
      m_density(density),                                      //
      m_viscosity(viscosity),                                  //
      m_airDrag(airDrag),                                      //
      m_internalRadiusMultiplier(1.),                          //
      m_ellipticalRadii(radiusA, radiusB),                     //
      m_baseRotation(baseRotation),                            //
      m_bendingMatrixBase(m_ellipticalRadii, m_baseRotation),  //
      m_youngsModulus(YoungsModulus),                          //
      m_shearModulus(shearModulus),                            //
      m_slipLength(0.0002),
      m_stretchMultiplier(stretchMultiplier),
      m_fixingMultiplier(fixingMultiplier),
      m_minBendingAngle(minBendingAngle),
      m_minBendingMultiplier(minBendingMultiplier),
      m_maxFlowGradientRatio(maxFlowGradientRatio),
      m_ks(m_ellipticalRadii, m_youngsModulus),  //
      m_kt(m_ellipticalRadii, m_shearModulus),
      m_flowOpenAtStart(flowOpenAtStart),
      m_flowOpenAtEnd(flowOpenAtEnd),
      m_collisionFree(false) {
  m_radiusMultiplier.resize(2);
  m_radiusMultiplier(0) = m_radiusMultiplier(1) = 1.;
  m_maxRadiusMultiplier = m_radiusMultiplier.maxCoeff();
}

ElasticStrandParameters::ElasticStrandParameters(
    int paramIndex, Scalar radiusA, Scalar radiusB, Scalar YoungsModulus,
    Scalar shearModulus, Scalar density, Scalar viscosity, Scalar airDrag,
    Scalar baseRotation, Scalar stretchMultiplier, Scalar fixingMultiplier,
    Scalar minBendingAngle, Scalar maxFlowGradientRatio,
    const VecXx& radiusMultiplier, Scalar minBendingMultiplier,
    bool flowOpenAtStart, bool flowOpenAtEnd)
    : m_paramIndex(paramIndex),
      m_density(density),                                      //
      m_viscosity(viscosity),                                  //
      m_airDrag(airDrag),                                      //
      m_internalRadiusMultiplier(1.),                          //
      m_ellipticalRadii(radiusA, radiusB),                     //
      m_baseRotation(baseRotation),                            //
      m_bendingMatrixBase(m_ellipticalRadii, m_baseRotation),  //
      m_youngsModulus(YoungsModulus),                          //
      m_shearModulus(shearModulus),                            //
      m_slipLength(0.0002),
      m_stretchMultiplier(stretchMultiplier),
      m_fixingMultiplier(fixingMultiplier),
      m_minBendingAngle(minBendingAngle),
      m_minBendingMultiplier(minBendingMultiplier),
      m_maxFlowGradientRatio(maxFlowGradientRatio),
      m_ks(m_ellipticalRadii, m_youngsModulus),  //
      m_kt(m_ellipticalRadii, m_shearModulus),
      m_radiusMultiplier(radiusMultiplier),
      m_flowOpenAtStart(flowOpenAtStart),
      m_flowOpenAtEnd(flowOpenAtEnd),
      m_collisionFree(false) {
  m_maxRadiusMultiplier = m_radiusMultiplier.maxCoeff();
}

ElasticStrandParameters& ElasticStrandParameters::operator=(
    const ElasticStrandParameters& other) {
  m_paramIndex = other.m_paramIndex;
  m_density = other.m_density;
  m_viscosity = other.m_viscosity;
  m_airDrag = other.m_airDrag;
  m_radiusMultiplier = other.m_radiusMultiplier;
  m_slipLength = other.m_slipLength;
  m_flowOpenAtStart = other.m_flowOpenAtStart;
  m_flowOpenAtEnd = other.m_flowOpenAtEnd;
  m_stretchMultiplier = other.m_stretchMultiplier;
  m_fixingMultiplier = other.m_fixingMultiplier;
  m_minBendingAngle = other.m_minBendingAngle;
  m_maxFlowGradientRatio = other.m_maxFlowGradientRatio;

  // Only the base dependencyNodes need to be copied
  m_ellipticalRadii.set(other.m_ellipticalRadii.get());
  m_baseRotation.set(other.m_baseRotation.get());
  m_youngsModulus.set(other.m_youngsModulus.get());
  m_shearModulus.set(other.m_shearModulus.get());

  m_maxRadiusMultiplier = m_radiusMultiplier.maxCoeff();

  return *this;
}

void ElasticStrandParameters::setBaseRotation(Scalar baseRotation) {
  m_baseRotation.set(baseRotation);
}

void ElasticStrandParameters::setDensity(Scalar density) {
  m_density = density;
}

void ElasticStrandParameters::setRadii(Scalar radiusA, Scalar radiusB) {
  m_ellipticalRadii.set(std::make_pair(radiusA, radiusB));
}

void ElasticStrandParameters::setShearModulus(Scalar shearModulus) {
  m_shearModulus.set(shearModulus);
}

void ElasticStrandParameters::setYoungsModulus(Scalar youngsModulus) {
  m_youngsModulus.set(youngsModulus);
}

void ElasticStrandParameters::setViscosity(Scalar viscosity) {
  m_viscosity = viscosity;
}

void ElasticStrandParameters::setAirDrag(Scalar airDrag) {
  m_airDrag = airDrag;
}

void ElasticStrandParameters::computeViscousForceCoefficients(Scalar dt) {
  const std::pair<Scalar, Scalar>& radii = m_ellipticalRadii.get();
  const Scalar m_radiusA = radii.first;
  const Scalar m_radiusB = radii.second;

  // Force coefficients are computed without the varying radius multiplier;
  // correct interpolation will be applied when they are accessed
  m_dt = dt;
  m_viscousKs = M_PI * m_radiusA * m_radiusB * 3 * m_viscosity / dt;
  m_viscousKt = M_PI_4 * m_radiusA * m_radiusB *
                (m_radiusA * m_radiusA + m_radiusB * m_radiusB) * m_viscosity /
                dt;
  m_viscousBendingCoefficientBase = 3 * m_viscosity / dt;
  m_viscousLamb = m_airDrag / dt;
}

Scalar ElasticStrandParameters::interpolatedRadiusMultiplier(
    int vtx, int num_verts) const {
  const Scalar s =
      vtx / (Scalar)(num_verts - 1) * (Scalar)(m_radiusMultiplier.size() - 1);
  const int i = (int)s;
  const Scalar frac = s - (Scalar)i;
  const int inext = std::min((int)m_radiusMultiplier.size() - 1, i + 1);

  return (frac * m_radiusMultiplier[i] +
          (1. - frac) * m_radiusMultiplier[inext]);
}

// Returns the multiplier e( theta ) such that the ellipse equation is
// r( theta ) = e( theta ) * majorRadius
// majorRadius being radiusA
Scalar ElasticStrandParameters::getRadiusShrinkAtAngle(
    const Scalar angle) const {
  const std::pair<Scalar, Scalar>& radii = m_ellipticalRadii.get();
  const Scalar m_radiusA = radii.first;
  const Scalar m_radiusB = radii.second;

  const Scalar a2 = m_radiusA * m_radiusA;
  const Scalar b2 = m_radiusB * m_radiusB;
  const Scalar ct = std::cos(angle);
  return m_radiusA / std::sqrt(b2 + (a2 - b2) * ct * ct);

  return radiusShrinkAtAngle(m_radiusA / m_radiusB, angle);
}

Scalar ElasticStrandParameters::radiusShrinkAtAngle(const Scalar ratio,
                                                    const Scalar angle) {
  const Scalar ct = std::cos(angle);
  return ratio / std::sqrt(1. + (ratio * ratio - 1.) * ct * ct);
}

}  // namespace strandsim
