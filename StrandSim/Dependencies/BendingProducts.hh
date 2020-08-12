/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BENDINGPRODUCTS_HH_
#define BENDINGPRODUCTS_HH_

#include "ElasticityParameters.hh"
#include "Kappas.hh"

namespace strandsim {

typedef std::vector<Mat2x, Eigen::aligned_allocator<Mat2x> >
    Mat2xArray;  ///< an array of 2d scalar matrices
typedef std::vector<Mat11x, Eigen::aligned_allocator<Mat11x> >
    Mat11xArray;  ///< an array of 11d scalar matrices

/**
 * \brief This class stores the products gradKappa^T B gradKappa that are used
 * in both the viscous and non-viscous bending forces.
 *
 * This product turned out to be one of the most costly operations in the
 * collision-free simulation.
 *
 * The matrix B is not the same in the viscous and non-viscous cases, but
 * differs only by a proportionality factor that can be applied when computing
 * the force Jacobian.
 *
 * Unit: cm^2
 */
class BendingProducts : public DependencyNode<Mat11xArray> {
 public:
  BendingProducts(BendingMatrixBase& bendingMatrixBase, GradKappas& gradKappas)
      : DependencyNode<Mat11xArray>(1, gradKappas.size()),  //
        m_bendingMatrixBase(bendingMatrixBase),             //
        m_gradKappas(gradKappas) {
#ifdef VERBOSE_DEPENDENCY_NODE
    std::cout << "Creating " << name() << ' ' << this << '\n';
#endif

    m_bendingMatrixBase.addDependent(this);
    m_gradKappas.addDependent(this);
  }

  virtual const char* name() const { return "BendingProducts"; }

 protected:
  virtual void compute();

  BendingMatrixBase& m_bendingMatrixBase;
  GradKappas& m_gradKappas;
};

}  // namespace strandsim

#endif /* BENDINGPRODUCTS_HH_ */
