/**
 * \copyright 2008 Robert Bridson
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef S_ARRAY_H
#define S_ARRAY_H

#include <algorithm>
#include <cassert>
#include <vector>

namespace strandsim {
namespace bridson {

template <class T>
struct Array {
  int nx, ny, nz;
  std::vector<T> a;

  Array(void) : nx(0), ny(0), nz(0) {}

  Array(int nx_, int ny_, int nz_)
      : nx(nx_), ny(ny_), nz(nz_), a(nx_ * ny_ * nz_) {}

  Array(int nx_, int ny_, int nz_, const T& value)
      : nx(nx_), ny(ny_), nz(nz_), a(nx_ * ny_ * nz_, value) {}

  void assign(int nx_, int ny_, int nz_, const T& value) {
    resize(nx_, ny_, nz_);
    assign(value);
  }

  void assign(const T& value) {
    for (unsigned int i = 0; i < a.size(); ++i) a[i] = value;
  }

  void clear(void) {
    nx = 0;
    ny = 0;
    nz = 0;
    a.clear();
  }

  bool empty(void) { return nx == 0 && ny == 0 && nz == 0; }

  void reserve(int nx_, int ny_, int nz_) { a.reserve(nx_ * ny_ * nz_); }

  void resize(int nx_, int ny_, int nz_) {
    nx = nx_;
    ny = ny_;
    nz = nz_;
    a.resize(nx * ny * nz);
  }

  void swap(Array& b) {
    std::swap(nx, b.nx);
    std::swap(ny, b.ny);
    std::swap(nz, b.nz);
    a.swap(b.a);
  }

  T* data(void) { return &a[0]; }

  T& operator()(int i, int j, int k) {
    assert(i < nx && j < ny && k < nz);
    return a[i + nx * j + nx * ny * k];
  }

  T& operator[](const int i) {
    assert(i >= 0 && i < nx * ny * nz);
    return a[i];
  }

  const T& operator()(int i, int j, int k) const {
    assert(i < nx && j < ny && k < nz);
    return a[i + nx * j + nx * ny * k];
  }

  const T& operator[](int i) const {
    assert(i >= 0 && i < nx * ny * nz);
    return a[i];
  }

  size_t size() const { return a.size(); }
};

typedef Array<double> Arrayd;
typedef Array<float> Arrayf;
typedef Array<long> Arrayl;
typedef Array<unsigned long> Arrayul;
typedef Array<int> Arrayi;
typedef Array<unsigned int> Arrayui;
typedef Array<short> Arrays;
typedef Array<unsigned short> Arrayus;
typedef Array<char> Arrayc;
typedef Array<unsigned char> Arrayuc;
// typedef Array<bool>           Arrayb; NOTE: this doesn't work as intended
// since std::vector<bool> is specialized

template <class T>
bool operator==(const Array<T>& a, const Array<T>& b) {
  return a.nx == b.nx && a.ny == b.ny && a.nz == b.nz && a.a == b.a;
}

template <class T>
bool operator!=(const Array<T>& a, const Array<T>& b) {
  return !(a == b);
}

template <class T>
bool operator<(const Array<T>& a, const Array<T>& b) {
  if (a.nx < b.nx)
    return true;
  else if (a.nx == b.nx) {
    if (a.ny < b.ny)
      return true;
    else if (a.ny == b.ny)
      if (a.nz < b.nz)
        return true;
      else if (a.nz == b.nz)
        return a.a < b.a;
  }
  return false;
}

template <class T>
bool operator>(const Array<T>& a, const Array<T>& b) {
  return b < a;
}

template <class T>
bool operator<=(const Array<T>& a, const Array<T>& b) {
  return !(b < a);
}

template <class T>
bool operator>=(const Array<T>& a, const Array<T>& b) {
  return !(a < b);
}

}  // namespace bridson
}  // namespace strandsim

#endif
