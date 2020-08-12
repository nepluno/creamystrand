/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef S_BOUNDINGBOX_HH_
#define S_BOUNDINGBOX_HH_

#include <stdint.h>

#include <algorithm>
#include <iostream>
#include <limits>

#include "Definitions.hh"

namespace strandsim {

typedef unsigned int uint;

template <typename ScalarT>
Eigen::Matrix<ScalarT, 3, 1> min(const Eigen::Matrix<ScalarT, 3, 1>& a,
                                 const Eigen::Matrix<ScalarT, 3, 1>& b) {
  return a.cwiseMin(b);
}

template <typename ScalarT>
Eigen::Matrix<ScalarT, 3, 1> max(const Eigen::Matrix<ScalarT, 3, 1>& a,
                                 const Eigen::Matrix<ScalarT, 3, 1>& b) {
  return a.cwiseMax(b);
}

template <typename ScalarT>
struct BoundingBox {
  typedef Eigen::Matrix<ScalarT, 3, 1> PointType;
  PointType min;
  PointType max;

  BoundingBox()
      : min(std::numeric_limits<ScalarT>::max(),
            std::numeric_limits<ScalarT>::max(),
            std::numeric_limits<ScalarT>::max()),
        max(-std::numeric_limits<ScalarT>::max(),
            -std::numeric_limits<ScalarT>::max(),
            -std::numeric_limits<ScalarT>::max()) {}

  explicit BoundingBox(PointType point) : min(point), max(point) {}

  BoundingBox(PointType minPoint, PointType maxPoint)
      : min(minPoint), max(maxPoint) {}

  BoundingBox(ScalarT minx, ScalarT miny, ScalarT minz, ScalarT maxx,
              ScalarT maxy, ScalarT maxz)
      : min(minx, miny, minz), max(maxx, maxy, maxz) {}

  virtual ~BoundingBox() {}

  bool isValid() const {
    return max[0] >= min[0] && max[1] >= min[1] && max[2] >= min[2];
  }

  void reset() { *this = BoundingBox(); }

  // Expand the box to contain the given point
  inline void insert(const PointType& p) {
    min = min.cwiseMin(p);
    max = max.cwiseMax(p);
  }

  // Expand the box to contain the given sphere
  inline void insert(const PointType& p, const ScalarT& radius) {
    min = min.cwiseMin(p - PointType::Constant(radius));
    max = max.cwiseMax(p + PointType::Constant(radius));
  }

  // Expand the box to contain the given box
  inline void insert(const BoundingBox& box) {
    min = min.cwiseMin(box.min);
    max = max.cwiseMax(box.max);
  }

  ScalarT volume() const { return (max - min).array().prod(); }

  ScalarT maxDim() const { return (max - min).maxCoeff(); }

  template <typename ScalarT1>
  friend std::ostream& operator<<(std::ostream& os,
                                  const BoundingBox<ScalarT1>& elem);

  template <typename ScalarT1>
  friend bool intersect(const BoundingBox<ScalarT1>& bbox_a,
                        const BoundingBox<ScalarT1>& bbox_b);

  template <typename ScalarT1>
  friend bool is_contained(const BoundingBox<ScalarT1>& bbox1,
                           const BoundingBox<ScalarT1>& bbox2,
                           const Scalar tol);

  template <typename ScalarT1>
  friend BoundingBox<ScalarT1> intersection(const BoundingBox<ScalarT1>& bbox1,
                                            const BoundingBox<ScalarT1>& bbox2);

  template <typename ScalarT1>
  friend bool is_left(const BoundingBox<ScalarT1>& bbox, const int axis,
                      const Scalar pivot);

  template <typename ScalarT1>
  friend void insert(BoundingBox<ScalarT1>& bbox,
                     const BoundingBox<ScalarT1>& bbox2);

  template <typename ScalarT1>
  friend BoundingBox<ScalarT1> merge(const BoundingBox<ScalarT1>& bbox1,
                                     const BoundingBox<ScalarT1>& bbox2);

  PointType ComputeCenter() const { return 0.5 * (min + max); }

 protected:
  virtual void print(std::ostream& os) const {
    os << "[" << min << " --- " << max << "]";
  }
};

template <typename ScalarT>
std::ostream& operator<<(std::ostream& os, const BoundingBox<ScalarT>& elem) {
  elem.print(os);
  return os;
}

template <typename ScalarT>
bool intersect(const BoundingBox<ScalarT>& bbox_a,
               const BoundingBox<ScalarT>& bbox_b) {
  if ((bbox_a.max[0] < bbox_b.min[0]) || (bbox_b.max[0] < bbox_a.min[0]) ||
      (bbox_a.max[1] < bbox_b.min[1]) || (bbox_b.max[1] < bbox_a.min[1]) ||
      (bbox_a.max[2] < bbox_b.min[2]) || (bbox_b.max[2] < bbox_a.min[2])) {
    return false;
  } else {
    return true;
  }
}

template <typename ScalarT>
bool is_contained(const BoundingBox<ScalarT>& bbox1,
                  const BoundingBox<ScalarT>& bbox2, const Scalar tol) {
  if (bbox1.max.x() > bbox2.max.x() + tol ||
      bbox1.min.x() < bbox2.min.x() - tol ||
      bbox1.max.y() > bbox2.max.y() + tol ||
      bbox1.min.y() < bbox2.min.y() - tol ||
      bbox1.max.z() > bbox2.max.z() + tol ||
      bbox1.min.z() < bbox2.min.z() - tol)
    return false;
  else
    return true;
}

template <typename ScalarT>
BoundingBox<ScalarT> intersection(const BoundingBox<ScalarT>& bbox1,
                                  const BoundingBox<ScalarT>& bbox2) {
  return BoundingBox<ScalarT>(max(bbox1.min, bbox2.min),
                              min(bbox1.max, bbox2.max));
}

template <typename ScalarT>
bool is_left(const BoundingBox<ScalarT>& bbox, const uint32_t axis,
             const ScalarT pivot) {
  return bbox.min[axis] + bbox.max[axis] < pivot * 2.f;
}

template <typename ScalarT>
void insert(BoundingBox<ScalarT>& bbox, const BoundingBox<ScalarT>& bbox2) {
  bbox.min = min(bbox.min, bbox2.min);
  bbox.max = max(bbox.max, bbox2.max);
}

template <typename ScalarT>
BoundingBox<ScalarT> merge(const BoundingBox<ScalarT>& bbox1,
                           const BoundingBox<ScalarT>& bbox2) {
  return BoundingBox<ScalarT>(min(bbox1.min, bbox2.min),
                              max(bbox1.max, bbox2.max));
}

}  // namespace strandsim

#endif /* S_BOUNDINGBOX_HH_ */
