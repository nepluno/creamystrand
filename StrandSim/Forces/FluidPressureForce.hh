/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSIM_FLUIDPRESSUREFORCE_HH
#define STRANDSIM_FLUIDPRESSUREFORCE_HH

#include "../Core/BandMatrixFwd.hh"
#include "ForceBase.hh"

namespace strandsim {

class FluidScriptingController;

class FluidPressureForce : public ForceBase {
 public:
  static const IndexType s_first =
      0;  // The first index on which this force can apply
  static const IndexType s_last = 0;  // The last index (counting from the end)

  typedef Vec3x LocalForceType;
  typedef Mat3x LocalJacobianType;
  typedef VecXx ForceVectorType;

  virtual std::string getName() const { return "fluid pressure"; }

  static void setScriptingController(
      const std::shared_ptr<FluidScriptingController>& controller);

  static void computeLocal(LocalForceType& localF, const ElasticStrand& strand,
                           StrandState& geometry, const IndexType vtx);

  static void computeLocal(LocalJacobianType& localJ,
                           const ElasticStrand& strand, StrandState& geometry,
                           const IndexType vtx);

  static void addInPosition(ForceVectorType& globalForce, const IndexType vtx,
                            const LocalForceType& localForce);

  static void addInPosition(JacobianMatrixType& globalJacobian,
                            const IndexType vtx,
                            const LocalJacobianType& localJacobian);

 protected:
  FluidPressureForce();
  virtual ~FluidPressureForce();

  static std::shared_ptr<FluidScriptingController> s_controller;
};

}  // namespace strandsim

#endif  // STRANDSIM_FLUIDDRAGFORCE_HH
