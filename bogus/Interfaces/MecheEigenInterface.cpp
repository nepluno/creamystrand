/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#define NOMINMAX
#ifdef WIN32
#include <Windows.h>
#endif

#include <algorithm>

#include "../../StrandSim/Core/BandMatrixFwd.hh"
#include "../../StrandSim/Core/Definitions.hh"
#include "../../StrandSim/Utils/SymmetricBandMatrixSolver.hh"
#include "../Core/Block.impl.hpp"
#include "../Core/Block.io.hpp"
#include "../Core/BlockSolvers/ADMM.hpp"
#include "../Core/BlockSolvers/Coloring.impl.hpp"
#include "../Core/BlockSolvers/GaussSeidel.hpp"
#include "../Core/BlockSolvers/ProductGaussSeidel.hpp"
#include "../Core/BlockSolvers/ProjectedGradient.hpp"
#include "../Core/Utils/Timer.hpp"
#include "../Extra/SecondOrder.impl.hpp"
#include "FrictionProblem.hpp"
#include "MecheEigenInterface.hpp"

#ifdef BOGUS_WITH_BOOST_SERIALIZATION
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <fstream>
#endif

namespace bogus {

MecheFrictionProblem::Options::Options()
    : maxThreads(0),
      maxIters(0),
      cadouxIters(0),
      tolerance(0),
      useInfinityNorm(false),
      algorithm(GaussSeidel),
      gsRegularization(0),
      gsColoring(false),
      gsSkipIters(-1),  // -1 means default
      tryZeroAsWell(true),
      ignoreVelocity(false),
      pgVariant(projected_gradient::SPG),
      admmProjStepSize(1),
      admmFpStepSize(1.e-3) {}

MecheFrictionProblem::MecheFrictionProblem()
    : m_primal(BOGUS_NULL_PTR(PrimalFrictionProblem<3u>)),
      m_dual(BOGUS_NULL_PTR(DualFrictionProblem<3u>)),
      m_lastSolveTime(0),
      m_f(BOGUS_NULL_PTR(double)),
      m_w(BOGUS_NULL_PTR(double)),
      m_mu(BOGUS_NULL_PTR(double)),
      m_out(&std::cout) {}

MecheFrictionProblem::~MecheFrictionProblem() { destroy(); }

void MecheFrictionProblem::destroy() {
  delete[] m_f;
  m_f = BOGUS_NULL_PTR(double);
  delete[] m_w;
  m_w = BOGUS_NULL_PTR(double);
  delete[] m_mu;
  m_mu = BOGUS_NULL_PTR(double);
  delete m_primal;
  m_primal = BOGUS_NULL_PTR(PrimalFrictionProblem<3u>);
  delete m_dual;
  m_dual = BOGUS_NULL_PTR(DualFrictionProblem<3u>);
}

void MecheFrictionProblem::ackCurrentResidual(unsigned GSIter, double err) {
#ifdef BOGUS_VERBOSE
  if (m_out) {
    *m_out << "Finished iteration " << GSIter << " with residual " << err
           << std::endl;
  }
#endif
  m_callback.trigger(GSIter, err, m_timer.elapsed());
}

void MecheFrictionProblem::reset() {
  destroy();

  m_primal = new PrimalFrictionProblem<3u>();
  m_lastSolveTime = 0;
}

void MecheFrictionProblem::fromPrimal(
    const int NObj,                    //!< number of subsystems
    const std::vector<unsigned> ndof,  //!< array of size \a NObj, the number of
                                       //!< degree of freedom of each subsystem
    const std::vector<strandsim::SymmetricBandMatrixSolver<double, 10>*>
        MassMat,  //!< the square ndof[i] long mass matrix of each subsystem
    const Eigen::VectorXd
        f_in,  //!< the constant term in \f$ M v + f= {}^t \! H (r - r_c) \f$
    const Eigen::VectorXd
        rc_in,  //!< the adhesion term in \f$ M v + f= {}^t \! H (r - r_c) \f$
    const Eigen::VectorXd filter_in,
    const int n_in,  //!< number of contact points
    const Eigen::VectorXd
        mu_in,  //!< array of size \a n giving the friction coeffs
    const Eigen::VectorXd
        yields_in,  //!< array of size \a n giving the friction coeffs
    const Eigen::VectorXd
        etas_in,  //!< array of size \a n giving the friction coeffs
    const Eigen::VectorXd
        powers_in,  //!< array of size \a n giving the friction coeffs
    const std::vector<Eigen::Matrix<double, 3, 3> >
        E_in,  // E matrix with form: 3*n_in x 3, !< array of size \f$ n \times
               // d \times d \f$ giving the \a n normals followed by the \a n
               // tangent vectors (and by again \a n tangent vectors if \a d is
               // 3). Said otherwise, \a E is a \f$ (nd) \times d \f$ matrix,
               // stored column-major, formed by \a n blocks of size \f$ d
               // \times d \f$ with each block being an orthogonal matrix (the
               // transition matrix from the world space coordinates \f$ (x_1,
               // x_2, x_3) \f$ to the local coordinates \f$ (x_N, x_{T1},
               // x_{T2}) \f$
    const Eigen::VectorXd w_in,  //!< array of size \a nd, the constant term in
                                 //!< \f$ u = H v + w \f$
    const int* const
        ObjA,  //!< array of size \a n, the first object involved in the \a i-th
               //!< contact (must be an internal object) (counted from 0)
    const int* const
        ObjB,  //!< array of size \a n, the second object involved in the \a
               //!< i-th contact (-1 for an external object) (counted from 0)
    const std::vector<strandsim::SparseRowMatx*>
        HA,  //!< array of size \a n, containing pointers to a dense,
             //!< colum-major matrix of size <c> d*ndof[ObjA[i]] </c>
             //!< corresponding to the H-matrix of <c> ObjA[i] </c>
    const std::vector<strandsim::SparseRowMatx*>
        HB,  //!< array of size \a n, containing pointers to a dense,
             //!< colum-major matrix of size <c> d*ndof[ObjA[i]] </c>
             //!< corresponding to the H-matrix of <c> ObjB[i] </c> (\c NULL for
             //!< an external object)
    const bool diagonalProblem) {
  reset();

  // Copy M
  // We don't actually need it after having computed a factorization of M, but
  // we keep it around in case we want to use dumpToFile()
  if (diagonalProblem) {
    m_primal->DiagM.reserve(NObj);
    m_primal->DiagM.setRows(ndof);
    m_primal->DiagM.setCols(ndof);
  } else {
    m_primal->M.reserve(NObj);
    m_primal->M.setRows(ndof);
    m_primal->M.setCols(ndof);
  }

  for (int i = 0; i < NObj; ++i) {
    if (diagonalProblem) {
      const int dim = MassMat[i]->matrix().rows();

      Eigen::DiagonalMatrix<double, Eigen::Dynamic> diagM(dim);

      for (int r = 0; r < dim; ++r) {
        diagM.diagonal()(r) = MassMat[i]->matrix()(r, r);
      }
      m_primal->DiagM.insertBack(i, i) = diagM;
    } else {
      Eigen::MatrixXd M(MassMat[i]->matrix().rows(),
                        MassMat[i]->matrix().cols());

      for (int r = 0; r < M.rows(); ++r) {
        for (int c = 0; c < M.cols(); ++c) {
          M(r, c) = MassMat[i]->matrix()(r, c);
        }
      }
      m_primal->M.insertBack(i, i) = M;
    }
  }

  if (diagonalProblem)
    m_primal->DiagM.finalize();
  else
    m_primal->M.finalize();

  // E
  m_primal->E.reserve(n_in);
  m_primal->E.setRows(n_in);
  m_primal->E.setCols(n_in);
  for (int i = 0; i < n_in; ++i) {
    m_primal->E.insertBack(i, i) = E_in[i];
  }
  m_primal->E.finalize();
  m_primal->E.cacheTranspose();

  // Build H
  m_primal->H.reserve(2 * n_in);
  m_primal->H.setRows(n_in);
  m_primal->H.setCols(ndof);

  for (std::ptrdiff_t i = 0; i < (std::ptrdiff_t)n_in; ++i) {
    const Eigen::Matrix3d Et = m_primal->E.diagonal(i).transpose();

    if (ObjB[i] == -1 && ObjA[i] != -1 && HA[i]) {
      m_primal->H.insert(i, ObjA[i]) = Et * Eigen::MatrixXd(*HA[i]);
    } else if (ObjA[i] == -1 && ObjB[i] != -1 && HB[i]) {
      m_primal->H.insert(i, ObjB[i]) = -Et * Eigen::MatrixXd(*HB[i]);
    } else if (ObjB[i] == ObjA[i] && ObjA[i] != -1 && HA[i] && HB[i]) {
      m_primal->H.insert(i, ObjA[i]) =
          Et * (Eigen::MatrixXd(*HA[i]) - Eigen::MatrixXd(*HB[i]));
    } else if (ObjA[i] != -1 && ObjB[i] != -1 && HA[i] && HB[i]) {
      m_primal->H.insert(i, ObjA[i]) = Et * Eigen::MatrixXd(*HA[i]);
      m_primal->H.insert(i, ObjB[i]) = -Et * Eigen::MatrixXd(*HB[i]);
    }
  }

  m_primal->H.finalize();

  m_primal->f = f_in;
  m_primal->rc = rc_in;
  m_primal->w = w_in;
  m_primal->mu = mu_in;
  m_primal->yields = yields_in;
  m_primal->etas = etas_in;
  m_primal->powers = powers_in;
  m_primal->filter = filter_in;

  if (diagonalProblem) {
    m_primal->computeDiagMInv();
  } else {
    m_primal->computeMInv();
  }
}

unsigned MecheFrictionProblem::nDegreesOfFreedom() const {
  return m_primal ? std::max(m_primal->M.rows(), m_primal->DiagM.rows()) : 0u;
}

unsigned MecheFrictionProblem::nContacts() const {
  return m_primal ? m_primal->H.rowsOfBlocks() : 0u;
}

void MecheFrictionProblem::computeDual(double regularization,
                                       bool diagonalProblem) {
  delete m_dual;
  m_dual = new DualFrictionProblem<3u>();
  m_dual->computeFrom(*m_primal, diagonalProblem);

  if (regularization > 0.) {
    m_dual->W += regularization * m_dual->W.Identity();
  }
}

double MecheFrictionProblem::solve(
    Eigen::VectorXd& r,  //!< length \a nd : initialization for \a r (in world
                         //!< space coordinates) + used to return computed r
    Eigen::VectorXd&
        v,  //!< length \a m: to return computed v ( or NULL if not needed )
    Eigen::VectorXd& r_world,
    int maxThreads,  //!< Maximum number of threads that the GS will use. If 0,
                     //!< use OpenMP default. If > 1, enable coloring to ensure
                     //!< deterministicity
    double tol,      //!< Gauss-Seidel tolerance. 0. means GS's default
    unsigned maxIters,   //!< Max number of iterations. 0 means GS's default
    bool staticProblem,  //!< If true, do not use DeSaxce change of variable
    bool herschelBulkleyProblem, bool doFrictionShrinking,
    double regularization,  //!< Coefficient to add to the diagonal of static
                            //!< problems / GS regularization coefficient for
                            //!< friction problems
    bool useInfinityNorm,  //!< Whether to use the infinity norm to evaluate the
                           //!< residual of the friction problem,
    bool useProjectedGradient,  //!< If true, use projected gradient algorithm
                                //!< instead of GS. Require either
                                //!< staticProblem=true, or cadouxIters > 0
    unsigned cadouxIters, bool diagonalProblem) {
  Options options;
  double problemRegularization = 0;

  options.maxThreads = maxThreads;
  options.maxIters = maxIters;
  options.cadouxIters = cadouxIters;

  options.tolerance = tol;
  options.useInfinityNorm = useInfinityNorm;

  if (useProjectedGradient) options.algorithm = ProjectedGradient;

  if (staticProblem) {
    problemRegularization = regularization;
  } else {
    options.gsRegularization = regularization;
  }

  return solve(r, v, r_world, options, staticProblem, herschelBulkleyProblem,
               doFrictionShrinking, problemRegularization, diagonalProblem);
}

double MecheFrictionProblem::solve(
    Eigen::VectorXd& r,  //!< length \a nd : initialization for \a r (in world
                         //!< space coordinates) + used to return computed r
    Eigen::VectorXd&
        v,  //!< length \a m: to return computed v ( or NULL if not needed )
    Eigen::VectorXd& r_world,
    const Options& options,  //!< Solver options
    bool staticProblem,  //!< If true, do not use DeSaxce change of variable, ie
                         //!< solve SOCQP -- useful for statics
    bool herschelBulkleyProblem, bool doFrictionShrinking,
    double problemRegularization, bool diagonalProblem) {
  assert(m_primal);
  const unsigned m = m_primal->H.cols();
  const unsigned n = m_primal->H.rowsOfBlocks();

  // r to local coords
  Eigen::VectorXd r_loc = m_primal->E.transpose() * r + m_primal->rc;
  ;

  double res;

  if (options.algorithm == ADMM || options.algorithm == DualAMA) {
    // Primal-dual solve
    m_timer.reset();

    if (options.algorithm == DualAMA) {
      bogus::DualAMA<bogus::PrimalFrictionProblem<3u>::HType> dama;
      dama.callback().connect(*this, &MecheFrictionProblem::ackCurrentResidual);

      if (options.tolerance != 0.) dama.setTol(options.tolerance);
      if (options.maxIters != 0) dama.setMaxIters(options.maxIters);
      dama.useInfinityNorm(options.useInfinityNorm);

      dama.setLineSearchIterations(0);
      dama.setFpStepSize(options.admmFpStepSize);
      dama.setProjStepSize(options.admmProjStepSize);

      res = m_primal->solveWith(dama, v.data(), r_loc.data(), staticProblem);
    } else {
      bogus::ADMM<bogus::PrimalFrictionProblem<3u>::HType> admm;
      admm.callback().connect(*this, &MecheFrictionProblem::ackCurrentResidual);

      admm.setStepSize(options.admmProjStepSize);

      if (options.tolerance != 0.) admm.setTol(options.tolerance);
      if (options.maxIters != 0) admm.setMaxIters(options.maxIters);
      admm.useInfinityNorm(options.useInfinityNorm);

      res = m_primal->solveWith(admm, 0., v.data(), r_loc.data());
    }

  } else {
    Signal<unsigned, double> callback;
    callback.connect(*this, &MecheFrictionProblem::ackCurrentResidual);

    if (options.algorithm == MatrixFreeGaussSeidel) {
      if (diagonalProblem) {
        typename PrimalFrictionProblem<3u>::DiagProductGaussSeidelType gs;
        if (options.tolerance != 0.) gs.setTol(options.tolerance);
        if (options.maxIters != 0) gs.setMaxIters(options.maxIters);
        if (options.gsSkipIters >= 0) gs.setSkipIters(options.gsSkipIters);

        gs.useInfinityNorm(options.useInfinityNorm);
        gs.setMaxThreads(options.maxThreads);
        gs.setAutoRegularization(options.gsRegularization);
        gs.doTryZeroAsWell(options.tryZeroAsWell);

        gs.callback().connect(callback);
        res = m_primal->solveWith(gs, r_loc.data(), staticProblem);
      } else {
        typename PrimalFrictionProblem<3u>::ProductGaussSeidelType gs;
        if (options.tolerance != 0.) gs.setTol(options.tolerance);
        if (options.maxIters != 0) gs.setMaxIters(options.maxIters);
        if (options.gsSkipIters >= 0) gs.setSkipIters(options.gsSkipIters);

        gs.useInfinityNorm(options.useInfinityNorm);
        gs.setMaxThreads(options.maxThreads);
        gs.setAutoRegularization(options.gsRegularization);
        gs.doTryZeroAsWell(options.tryZeroAsWell);

        gs.callback().connect(callback);
        res = m_primal->solveWith(gs, r_loc.data(), staticProblem);
      }
    } else {
      // If dual has not been computed yet
      if (!m_dual) {
        computeDual(problemRegularization, diagonalProblem);
      }

      // Proper solving

      m_timer.reset();
      if (options.algorithm == ProjectedGradient) {
        DualFrictionProblem<3u>::ProjectedGradientType pg;
        if (options.tolerance != 0.) pg.setTol(options.tolerance);
        if (options.maxIters != 0) pg.setMaxIters(options.maxIters);

        pg.useInfinityNorm(options.useInfinityNorm);
        pg.setDefaultVariant(options.pgVariant);

        if (staticProblem || options.cadouxIters == 0) {
          pg.callback().connect(callback);
          res = m_dual->solveWith(pg, r_loc.data(), staticProblem,
                                  herschelBulkleyProblem);
        } else {
          res = m_dual->solveCadoux(pg, r_loc.data(), options.cadouxIters,
                                    &callback);
        }
      } else {
        // Setup GS parameters
        bogus::DualFrictionProblem<3u>::GaussSeidelType gs;

        if (options.tolerance != 0.) gs.setTol(options.tolerance);
        if (options.maxIters != 0) gs.setMaxIters(options.maxIters);
        if (options.gsSkipIters >= 0) gs.setSkipIters(options.gsSkipIters);

        gs.setMaxThreads(options.maxThreads);
        gs.setAutoRegularization(options.gsRegularization);
        gs.useInfinityNorm(options.useInfinityNorm);

        m_dual->undoPermutation();

        const bool useColoring = options.maxThreads != 1 && options.gsColoring;
        gs.coloring().update(useColoring, m_dual->W);

        if (useColoring) {
          m_dual->applyPermutation(gs.coloring().permutation);
          gs.coloring().resetPermutation();
        }
        m_dual->W.cacheTranspose();

        if (staticProblem || options.cadouxIters == 0) {
          gs.doTryZeroAsWell(options.tryZeroAsWell);
          gs.callback().connect(callback);
          res = m_dual->solveWith(gs, r_loc.data(), staticProblem,
                                  herschelBulkleyProblem);
        } else {
          gs.doTryZeroAsWell(false);
          res = m_dual->solveCadoux(gs, r_loc.data(), options.cadouxIters,
                                    &callback);
        }
      }
    }
  }

  m_lastSolveTime = m_timer.elapsed();
#ifdef BOGUS_VERBOSE
  if (m_out && n != 0) {
    *m_out << "Max coeff: " << r_loc.lpNorm<Eigen::Infinity>() << std::endl;
  }
#endif
  Eigen::VectorXd r_new = r_loc - m_primal->rc;

  if (doFrictionShrinking) {
    for (unsigned int i = 0; i < n; ++i) {
      const double nrt = r_new.segment<2>(i * 3 + 1).norm();
      if (nrt < 1e-20) continue;

      const double prop =
          std::max(0., nrt - m_primal->mu[i] * m_primal->rc(i * 3)) / nrt;
      r_new.segment<2>(i * 3 + 1) *= prop;
    }
  }

  if (herschelBulkleyProblem) {
    for (unsigned int i = 0; i < n; ++i) {
      r_new(i * 3) = 0.0;
    }
  }

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    const double maxImpulseNorm = m_primal->filter[i];

    if (maxImpulseNorm > 0.) {
      const double len = r_new.segment<3>(i * 3).norm();
      if (len > maxImpulseNorm) {
        r_new.segment<3>(i * 3) =
            r_new.segment<3>(i * 3).normalized() * maxImpulseNorm;
      }
    }
  }

  r_world = m_primal->H.transpose() * r_new;

  if (!options.ignoreVelocity) {
    if (diagonalProblem)
      v = m_primal->DiagMInv * (r_world - m_primal->f);
    else
      v = m_primal->MInv * (r_world - m_primal->f);
  }

  // r to world coords
  r = m_primal->E * (r_loc - m_primal->rc);

  return res;
}

void MecheFrictionProblem::setOutStream(std::ostream* out) { m_out = out; }
}  // namespace bogus
