/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "BendingProducts.hh"

#include "../Core/ElasticStrandUtils.hh"

namespace strandsim {

void BendingProducts::compute() {
  m_value.resize(m_size);
  const Mat2x& bendingMatrix = m_bendingMatrixBase.get();
  const GradKArrayType& gradKappas = m_gradKappas.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    symBProduct<11>(m_value[vtx], bendingMatrix,
                    gradKappas[vtx].block<11, 2>(0, 0));
    symBProductAdd<11>(m_value[vtx], bendingMatrix,
                       gradKappas[vtx].block<11, 2>(0, 2));
  }

  setDependentsDirty();
}

}  // namespace strandsim
