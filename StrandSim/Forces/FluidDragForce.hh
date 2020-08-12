/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSIM_FLUIDDRAGFORCE_HH
#define STRANDSIM_FLUIDDRAGFORCE_HH

#include "../Core/BandMatrixFwd.hh"
#include "ForceBase.hh"

namespace strandsim {

class FluidScriptingController;

class FluidDragForce : public ForceBase {
 public:
  static const IndexType s_first =
      0;  // The first index on which this force can apply
  static const IndexType s_last = 0;  // The last index (counting from the end)

  typedef Vec3x LocalForceType;
  typedef Mat3x LocalJacobianType;
  typedef VecXx ForceVectorType;

  virtual std::string getName() const { return "fluid drag"; }

  static void setScriptingController(
      const std::shared_ptr<FluidScriptingController>& controller);

  static void computeLocal(LocalForceType& localF, const ElasticStrand& strand,
                           StrandState& geometry, const IndexType vtx);

  static void computeLocal(LocalJacobianType& localJ,
                           const ElasticStrand& strand, StrandState& geometry,
                           const IndexType vtx);

  static Scalar computeLocal(LocalForceType& localF, const Vec3x& edge,
                             const Vec3x& u_s0, const Vec3x& u_s,
                             const Vec3x& u_f, const Scalar intersectLength,
                             const Scalar radius, const Scalar rho,
                             const Scalar viscosity, const Scalar behavior,
                             const Scalar yield, const Scalar epsilon,
                             bool useConstantDrag);

  static void computeLocal(LocalJacobianType& localJ, const Vec3x& edge,
                           const Vec3x& u_s0, const Vec3x& u_f,
                           const Scalar intersectLength, const Scalar radius,
                           const Scalar rho, const Scalar viscosity,
                           const Scalar behavior, const Scalar yield,
                           const Scalar epsilon, const Scalar dt,
                           bool useConstantDrag);

  static Scalar computeCoeff(const Vec3x& edge, const Vec3x& u_s0,
                             const Vec3x& u_f0, const Scalar intersectLength,
                             const Scalar radius, const Scalar rho,
                             const Scalar viscosity, const Scalar behavior,
                             const Scalar yield, const Scalar epsilon,
                             bool useConstantDrag);

  static void addInPosition(ForceVectorType& globalForce, const IndexType vtx,
                            const LocalForceType& localForce);

  static void addInPosition(JacobianMatrixType& globalJacobian,
                            const IndexType vtx,
                            const LocalJacobianType& localJacobian);

 protected:
  FluidDragForce();
  virtual ~FluidDragForce();

  static std::shared_ptr<FluidScriptingController> s_controller;
};

}  // namespace strandsim

#endif  // STRANDSIM_FLUIDDRAGFORCE_HH
