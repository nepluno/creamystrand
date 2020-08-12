/**
 * \copyright 2015 Xinxin Zhang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "GeometricLevelGen.hh"

//#define AMG_VERBOSE

namespace strandsim {

template <class T>
void levelGen<T>::generateRP(const FixedSparseMatrix<T> &A,
                             FixedSparseMatrix<T> &R, FixedSparseMatrix<T> &P,
                             int ni, int nj, int nk) {
  int nni = std::ceil((double)ni / 2.0);
  int nnj = std::ceil((double)nj / 2.0);
  int nnk = std::ceil((double)nk / 2.0);
  SparseMatrix<T> r;
  SparseMatrix<T> p;
  p.resize(ni * nj * nk);
  p.zero();
  r.resize(nni * nnj * nnk);
  r.zero();

  for (int k = 0; k < nnk; k++)
    for (int j = 0; j < nnj; j++)
      for (int i = 0; i < nni; i++) {
        unsigned int index = (k * nnj + j) * nni + i;
        for (int kk = 0; kk <= 1; kk++)
          for (int jj = 0; jj <= 1; jj++)
            for (int ii = 0; ii <= 1; ii++) {
              int iii = i * 2 + ii;
              int jjj = j * 2 + jj;
              int kkk = k * 2 + kk;
              if (iii < ni && jjj < nj && kkk < nk) {
                unsigned int index2 = (kkk * nj + jjj) * ni + iii;
                r.set_element(index, index2, (T)0.125);
                p.set_element(index2, index, 1.0);
              }
            }
      }

  R.construct_from_matrix(r);
  P.construct_from_matrix(p);
  r.clear();
  p.clear();

  // transposeMat(R,P,(T)8.0);
}

template <class T>
void levelGen<T>::generateLevelsGalerkinCoarsening(
    vector<FixedSparseMatrix<T> > &A_L, vector<FixedSparseMatrix<T> > &R_L,
    vector<FixedSparseMatrix<T> > &P_L, vector<Vec3i> &S_L, int &total_level,
    FixedSparseMatrix<T> &A, int ni, int nj, int nk) {
#ifdef AMG_VERBOSE
  cout << "building levels ...... " << endl;
#endif
  A_L.resize(0);
  R_L.resize(0);
  P_L.resize(0);
  S_L.resize(0);
  total_level = 1;
  A_L.push_back(std::move(A));
  S_L.push_back(Vec3i(ni, nj, nk));
  int nni = ni, nnj = nj, nnk = nk;
  unsigned int unknowns = ni * nj * nk;
  while (unknowns > 16 * 16 * 16) {
    A_L.push_back(FixedSparseMatrix<T>());
    R_L.push_back(FixedSparseMatrix<T>());
    P_L.push_back(FixedSparseMatrix<T>());
    nni = std::ceil((double)nni / 2.0);
    nnj = std::ceil((double)nnj / 2.0);
    nnk = std::ceil((double)nnk / 2.0);

    S_L.push_back(Vec3i(nni, nnj, nnk));
    unknowns = nni * nnj * nnk;
    total_level++;
  }

  for (int i = 0; i < total_level - 1; i++) {
    generateRP((A_L[i]), (R_L[i]), (P_L[i]), S_L[i][0], S_L[i][1], S_L[i][2]);
    FixedSparseMatrix<T> temp;
    multiplyMat((A_L[i]), (P_L[i]), temp, 1.0);
    multiplyMat((R_L[i]), temp, (A_L[i + 1]), 0.5);
    temp.resize(0);
    temp.clear();
  }
#ifdef AMG_VERBOSE
  cout << "build levels done" << endl;
#endif
}

template struct levelGen<double>;
};  // namespace strandsim
