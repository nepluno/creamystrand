/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_COLORING_HPP
#define BOGUS_COLORING_HPP

#include "../Block.fwd.hpp"

namespace bogus {

//! Coloring algorithm to determine which rows of a matrix can be treated in
//! parallel
/*! Computes a permutation of the rows indices so that they become contiguous
 * for each color */
struct Coloring {
  //! Computed permuation so that each color is contiguous
  std::vector<std::size_t> permutation;
  //! Index of first row for each color
  std::vector<std::ptrdiff_t> colors;

  Coloring() {}

  //! Computes a coloring for \p matrix, or simply reset it if \p enable is
  //! false
  template <typename Derived>
  void update(const bool enable, const BlockMatrixBase<Derived>& matrix);

  std::size_t size() const { return permutation.size(); }

  //! Sets the permutation to the identity. Keep the current colors.
  void resetPermutation() {
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
    for (std::ptrdiff_t i = 0; i < (std::ptrdiff_t)permutation.size(); ++i) {
      permutation[i] = i;
    }
  }

 private:
  void reset(std::size_t n) {
    colors.clear();
    colors.push_back(0);
    colors.push_back(n);

    permutation.resize(n);
    resetPermutation();
  }

  template <typename Derived>
  void compute(const SparseBlockMatrixBase<Derived>& matrix);

  template <typename Derived>
  void compute(const BlockMatrixBase<Derived>& matrix);
};

}  // namespace bogus

#endif
