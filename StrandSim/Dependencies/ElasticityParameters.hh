/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ELASTICPARAMETERS_HH_
#define ELASTICPARAMETERS_HH_

#include <math.h>
#ifndef M_PI
#define M_PI 3.141592653589793238462650288
#endif
#ifndef M_PI_4
#define M_PI_4 M_PI / 4
#endif
#include "DependencyNode.hh"

namespace strandsim {

/**
 * Unit: cm
 */
class EllipticalRadii : public DependencyNode<std::pair<Scalar, Scalar> > {
 public:
  EllipticalRadii(Scalar radiusA, Scalar radiusB)
      : DependencyNode<std::pair<Scalar, Scalar> >(
            std::make_pair(radiusA, radiusB)) {
#ifdef VERBOSE_DEPENDENCY_NODE
    std::cout << "Creating " << name() << ' ' << this << '\n';
#endif

    setClean();
  }

  virtual const char* name() const { return "EllipticalRadii"; }

  friend class boost::serialization::access;
  template <class Archive>
  void save(Archive& ar, const int version) const {
    // Can't use get() or compute() here as boost require save() to be const.
    // So we'll just issue a warning if someone tries to save a dirty value.
    if (isDirty()) {
      WarningStream(g_log, "") << "Saving dirty value for " << name();
    }
    ar& m_value.first;
    ar& m_value.second;
  }
  template <class Archive>
  void load(Archive& ar, const int version) {
    std::pair<Scalar, Scalar> value;
    ar& value.first;
    ar& value.second;
    set(value);
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()

 protected:
  virtual void compute() {}
};

/**
 * Unit: no dimension
 */
class BaseRotation : public DependencyNode<Scalar> {
 public:
  BaseRotation(Scalar baseRotation) : DependencyNode<Scalar>(baseRotation) {
#ifdef VERBOSE_DEPENDENCY_NODE
    std::cout << "Creating " << name() << ' ' << this << '\n';
#endif

    setClean();
  }

  virtual const char* name() const { return "BaseRotation"; }

 protected:
  virtual void compute() {}
};

/**
 * \brief This contains the bending matrix base, must be multiplied by the
 * appropriate viscous or non-viscous coefficient (with optional interpolation
 * factor).
 *
 * Unit: cm^4
 */
class BendingMatrixBase : public DependencyNode<Mat2x> {
 public:
  BendingMatrixBase(EllipticalRadii& ellipticalRadii,
                    BaseRotation& baseRotation)
      : DependencyNode<Mat2x>(Mat2x()),      //
        m_ellipticalRadii(ellipticalRadii),  //
        m_baseRotation(baseRotation) {
#ifdef VERBOSE_DEPENDENCY_NODE
    std::cout << "Creating " << name() << ' ' << this << '\n';
#endif

    m_ellipticalRadii.addDependent(this);
    m_baseRotation.addDependent(this);
  }

  virtual const char* name() const { return "BendingMatrixBase"; }

 protected:
  virtual void compute() {
    const std::pair<Scalar, Scalar>& radii = m_ellipticalRadii.get();
    const Scalar baseRotation = m_baseRotation.get();

    Mat2x& B = m_value;
    B(0, 0) = M_PI_4 * radii.second * cube(radii.first);
    B(1, 1) = M_PI_4 * radii.first * cube(radii.second);
    // rotate cross section by a constant angle
    const Mat2x& rot =
        Eigen::Rotation2D<Scalar>(baseRotation).toRotationMatrix();
    B = rot * B * rot.transpose();
    B(0, 1) = B(1, 0) =
        0.5 * (B(0, 1) + B(1, 0));  // For perfect numerical symmetry

    setDependentsDirty();
  }

  EllipticalRadii& m_ellipticalRadii;
  BaseRotation& m_baseRotation;
};

/**
 * Unit: dPa = g cm^-1 s^-2
 */
class YoungsModulus : public DependencyNode<Scalar> {
 public:
  YoungsModulus(Scalar youngsModulus) : DependencyNode<Scalar>(youngsModulus) {
    setClean();
  }

  virtual const char* name() const { return "YoungsModulus"; }

 protected:
  virtual void compute() {}
};

/**
 * Unit: dPa = g cm^-1 s^-2
 */
class ShearModulus : public DependencyNode<Scalar> {
 public:
  ShearModulus(Scalar shearModulus) : DependencyNode<Scalar>(shearModulus) {
    setClean();
  }

  virtual const char* name() const { return "ShearModulus"; }

 protected:
  virtual void compute() {}
};

/**
 * Unit: 10^-5 N = g cm s^-2
 */
class ElasticKs : public DependencyNode<Scalar> {
 public:
  ElasticKs(EllipticalRadii& ellrad, YoungsModulus& ym)
      : DependencyNode<Scalar>(
            std::numeric_limits<Scalar>::signaling_NaN()),  //
        m_ellipticalRadii(ellrad),                          //
        m_youngsModulus(ym)                                 //
  {
    m_ellipticalRadii.addDependent(this);
    m_youngsModulus.addDependent(this);
  }

  virtual const char* name() const { return "ElasticKs"; }

 protected:
  virtual void compute() {
    const std::pair<Scalar, Scalar>& radii = m_ellipticalRadii.get();
    const Scalar youngsModulus = m_youngsModulus.get();

    m_value = M_PI * radii.first * radii.second * youngsModulus;

    setDependentsDirty();
  }

  EllipticalRadii& m_ellipticalRadii;
  YoungsModulus& m_youngsModulus;
};

/**
 * Unit: 10^-5 cm^2 N = g cm^3 s^-2
 */
class ElasticKt : public DependencyNode<Scalar> {
 public:
  ElasticKt(EllipticalRadii& ellrad, ShearModulus& sm)
      : DependencyNode<Scalar>(
            std::numeric_limits<Scalar>::signaling_NaN()),  //
        m_ellipticalRadii(ellrad),                          //
        m_shearModulus(sm) {
    m_ellipticalRadii.addDependent(this);
    m_shearModulus.addDependent(this);
  }

  virtual const char* name() const { return "ElasticKt"; }

 protected:
  virtual void compute() {
    const std::pair<Scalar, Scalar>& radii = m_ellipticalRadii.get();
    const Scalar shearModulus = m_shearModulus.get();

    m_value = M_PI_4 * radii.first * radii.second *
              (radii.first * radii.first + radii.second * radii.second) *
              shearModulus;

    setDependentsDirty();
  }

  EllipticalRadii& m_ellipticalRadii;
  ShearModulus& m_shearModulus;
};

}  // namespace strandsim

#endif /* ELASTICPARAMETERS_HH_ */
