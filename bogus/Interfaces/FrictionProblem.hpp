/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * So-bogus is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * So-bogus is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with So-bogus.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef BOGUS_FRICTION_PROBLEM_HPP
#define BOGUS_FRICTION_PROBLEM_HPP

#include "../Core/Block.hpp"
#include "../Core/BlockSolvers.fwd.hpp"
#include "../Core/Utils/Signal.hpp"
#include "../Extra/SecondOrder.fwd.hpp"

namespace bogus {

//! Primal Coulomb friction problem for a block-diagonal mass matrix M with
//! dense blocks
/*!
           \f[
                \left\{
                  \begin{array}{rcl}
                        M v &=& H^T r - f \\
                        u &=& H v + E^T w \\
                        &s.t.& law (x,y)
                  \end{array}
                \right.
          \f]
        where law is either SOCLaw (meaning solving a SOCQP) or CoulombLaw,
        and \p v are the velocities, \p r the contact forces and \p u the
   relative velocities at contact points. E is a rotation matrix transforming
   local contact basis coordinates into world coordinates , is is only useful
   for consistency with legacy interface.

        Useful for solving friction problems in discrete systems with reduced
   coordinate models. Solving may be performed directly on this formulation, or
   on the dual formulation by constructing a DualFrictionProblem from this
   object. See also \ref block_solvers_ns .
*/
template <unsigned Dimension>
struct PrimalFrictionProblem {
  // Primal Data

  typedef SparseBlockMatrix<Eigen::MatrixXd> MType;
  //! M -- mass matrix
  MType M;

  typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagonalXd;
  typedef SparseBlockMatrix<DiagonalXd> DiagMType;

  DiagMType DiagM;

  //! E -- local rotation matrix ( contact basis coordinates to world
  //! coordinates )
  bogus::SparseBlockMatrix<Eigen::Matrix<double, Dimension, Dimension> > E;

  typedef Eigen::Matrix<double, Dimension, Eigen::Dynamic> HBlock;
  typedef SparseBlockMatrix<HBlock, UNCOMPRESSED> HType;
  //! H -- deformation gradient \f$ \frac{\partial u}{\partial v} \f$ (
  //! generalized coordinates <-> contact basis coordinates )
  HType H;

  //! Adhesion impulses ( with unit force * time, in generalized coordinate )
  Eigen::VectorXd rc;
  //! External forces
  Eigen::VectorXd f;
  //! Free velocity ( such that u = Hv + w )
  Eigen::VectorXd w;
  //! Coulomb friction coefficients
  Eigen::VectorXd mu;

  Eigen::VectorXd yields;

  Eigen::VectorXd etas;

  Eigen::VectorXd powers;

  Eigen::VectorXd filter;

  // Cached data

  //! M^-1
  typedef SparseBlockMatrix<LU<Eigen::MatrixBase<Eigen::MatrixXd> > > MInvType;
  MInvType MInv;

  typedef SparseBlockMatrix<DiagonalXd> DiagMInvType;
  DiagMInvType DiagMInv;

  //! Computes MInv from M. Required to build a DualFrictionProblem for the
  //! PrimalFrictionProblem, or to use the ADMM and matrix-free Gauss-Seidel
  //! solvers.
  void computeMInv();

  void computeDiagMInv();
  // Primal-dual solve functions

  typedef ADMM<HType> ADMMType;
  typedef DualAMA<HType> DualAMAType;
  typedef ProductGaussSeidel<HType, MInvType, true> ProductGaussSeidelType;
  typedef ProductGaussSeidel<HType, DiagMInvType, true>
      DiagProductGaussSeidelType;

  //! Matrix-free Gauss-Seidel solver
  /*! \note Requires the computation of a factorization of M with computeMInv()
    \param r Resulting forces and initial guess (in contact basis coordinates)
    \param staticProblem If true, solve this problem as a \b SOCQP instead of a
    Coulomb Friction problem
  */
  double solveWith(ProductGaussSeidelType &pgs, double *r,
                   const bool staticProblem = false) const;

  double solveWith(DiagProductGaussSeidelType &pgs, double *r,
                   const bool staticProblem = false) const;
  //! ADMM (Alternating Direction Method of Multipliers) on the primal objective
  //! function.
  /*!
    May only be used to solve SOCQP (static problems), not Coulomb friction
    problems \param lambda proximal coefficient of the quadratic part (lambda =
    0 means AMA) \param v Resulting velocities and initial guess \param r
    Resulting forces and initial guess (in contact basis coordinates) \note
    Requires the computation of a factorization of M with computeMInv()
  */
  double solveWith(ADMMType &admm, double lambda, double *v, double *r) const;
  //! AMA (Alternating Minimization Algorithm) on the dual objective function. .
  /*! Does not require a factorization of M
    \param v Resulting velocities and initial guess
    \param r Resulting forces and initial guess (in contact basis coordinates)
    \param staticProblem If true, solve this problem as a \b SOCQP instead of a
    Coulomb Friction problem
  */
  double solveWith(DualAMAType &dama, double *v, double *r,
                   const bool staticProblem = false) const;
};

//! Dual representation of a Coulomb friction problem, with explicit Delassus
//! operator
/*!
  u = W r + b such that  u and r satisfy CoulombLaw(mu) or SOCLaw(mu),
  where \p W is a symmetric, positive semi-definite matrix with \f$ d \times d
  \f$ blocks, \p u are the relative velocities and \p r are the contact forces.
  May be constructed from a PrimalFrictionProblem, or by directly initializing
  each data member to accomodate problems with different mass matrix structures.
  See also \ref block_solvers_ns .
*/
template <unsigned Dimension>
struct DualFrictionProblem {
  typedef SparseBlockMatrix<
      Eigen::Matrix<double, Dimension, Dimension, Eigen::RowMajor>, SYMMETRIC>
      WType;

  typedef GaussSeidel<WType> GaussSeidelType;
  typedef ProjectedGradient<WType> ProjectedGradientType;

  typedef SOCLaw<Dimension, double, true> CoulombLawType;
  typedef SOCLaw<Dimension, double, false> SOCLawType;

  typedef HBLaw<Dimension, double> HBLawType;

  typedef Signal<unsigned, double> SignalType;

  //! W -- Delassus operator
  WType W;

  //! Rhs ( such that u = Wr + b )
  Eigen::VectorXd b;

  //! Coulomb friction coefficients
  Eigen::VectorXd mu;

  Eigen::VectorXd yields;

  Eigen::VectorXd etas;

  Eigen::VectorXd powers;

  //! Computes this DualFrictionProblem from the given \p primal
  /*! \warning Assumes MInv has been computed */
  void computeFrom(const PrimalFrictionProblem<Dimension> &primal,
                   const bool diagonalProblem = false);

  //! Solves this problem
  /*!
    \param gs The GaussSeidel< WType > solver to use
    \param r  Both the initial guess and the result
    \param staticProblem If true, solve this problem as a \b SOCQP instead of a
    Coulomb Friction problem \returns the error as returned by the
    GaussSeidel::solve() function
    */
  double solveWith(GaussSeidelType &gs, double *r,
                   const bool staticProblem = false,
                   const bool herschelBulkleyProblem = false) const;
  //! Same as above
  /*! \warning staticProblem defaults tp true (as solving Coulomb probles with
   * PG is unreliable)*/
  double solveWith(ProjectedGradientType &pg, double *r,
                   const bool staticProblem = true,
                   const bool herschelBulkleyProblem = false) const;

  //! Evaluate a residual using the GS's error function
  /*!
    \param gs The GaussSeidel< WType > solver to use
    \param r  Both the current force
    \param staticProblem If true, eval this problem as a \b SOCQP instead of a
    Coulomb Friction problem

    \returns the error as returned by the GaussSeidel::eval() function
    */
  double evalWith(const GaussSeidelType &gs, const double *r,
                  const bool staticProblem = false,
                  const bool herschelBulkleyProblem = false) const;
  //! Same as above
  /*! \warning staticProblem defaults tp true (as solving Coulomb probles with
   * PG is unreliable)*/
  double evalWith(const ProjectedGradientType &gs, const double *r,
                  const bool staticProblem = true,
                  const bool herschelBulkleyProblem = false) const;

  //! Solves this problem using the Cadoux algorithm ( with fixed-point
  //! iteration )
  /*!
    See \cite ACLM11
    \param gs The GaussSeidel< WType > solver to use
    \param r  Both the initial guess and the result
    \param fpIterations Number of fixed-point iterations
    \param callback 0, or a pointer to a user-defined function that takes (
    unsigned iteration, double residual ) as arguments \returns the error as
    returned by the GaussSeidel::solve() function
    */
  double solveCadoux(
      GaussSeidelType &gs, double *r, const unsigned fpIterations,
      const SignalType *callback = BOGUS_NULL_PTR(const SignalType)) const;
  double solveCadoux(
      ProjectedGradientType &pg, double *r, const unsigned fpIterations,
      const SignalType *callback = BOGUS_NULL_PTR(const SignalType)) const;

  //! Idem as solveCadoux, but interpreting the problem as r = Wu + b
  double solveCadouxVel(
      GaussSeidelType &gs, double *u, const unsigned fpIterations,
      const SignalType *callback = BOGUS_NULL_PTR(const SignalType)) const;
  double solveCadouxVel(
      ProjectedGradientType &pg, double *u, const unsigned fpIterations,
      const SignalType *callback = BOGUS_NULL_PTR(const SignalType)) const;

  //! Apply a permutation on the contact indices
  /*!
   Useful for achieving better memory locality when using the Coloring
   functionality of the GaussSeidel algorithm. \sa Coloring
   * \warning To use the permutation releated functions, all the blocks have to
   have the same size
  */
  void applyPermutation(const std::vector<std::size_t> &permutation);
  void undoPermutation();
  bool permuted() const { return !m_permutation.empty(); }

  const std::vector<std::size_t> &permutation() const { return m_permutation; }
  const std::vector<std::size_t> &invPermutation() const {
    return m_invPermutation;
  }

 private:
  // Current permutation of contact indices
  std::vector<std::size_t> m_permutation;
  std::vector<std::size_t> m_invPermutation;
};

}  // namespace bogus

#endif
