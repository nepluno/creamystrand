/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "MathUtilities.hh"

#include <iomanip>

#include "../Core/Definitions.hh"
#include "ThreadUtils.hh"

namespace strandsim {
bool approxSymmetric(const MatXx &A, const Scalar &eps) {
  for (int i = 0; i < A.rows(); ++i)
    for (int j = i + 1; j < A.cols(); ++j)
      if (fabs(A(i, j) - A(j, i)) >= eps) return false;
  return true;
}

Scalar ScalarRand(const Scalar min, const Scalar max) {
  static thread_local std::mt19937 generator;
  std::uniform_real_distribution<Scalar> distribution(min, max);
  return distribution(generator);
}

Scalar point_triangle_distance(const Vec3x &p, const Vec3x &a, const Vec3x &b,
                               const Vec3x &c, Scalar &t1, Scalar &t2,
                               Scalar &t3) {
  Scalar ab[3], ac[3], ap[3], bp[3];

  ab[0] = b[0] - a[0];
  ab[1] = b[1] - a[1];
  ab[2] = b[2] - a[2];

  ac[0] = c[0] - a[0];
  ac[1] = c[1] - a[1];
  ac[2] = c[2] - a[2];

  ap[0] = p[0] - a[0];
  ap[1] = p[1] - a[1];
  ap[2] = p[2] - a[2];

  Scalar d1 = ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2];
  Scalar d2 = ac[0] * ap[0] + ac[1] * ap[1] + ac[2] * ap[2];

  if ((d1 <= 0.0f) && (d2 <= 0.0f)) {
    t1 = 1.0f;
    t2 = 0.0f;
    t3 = 0.0f;

    return std::sqrt((p[0] - a[0]) * (p[0] - a[0]) +
                     (p[1] - a[1]) * (p[1] - a[1]) +
                     (p[2] - a[2]) * (p[2] - a[2]));
  }

  bp[0] = p[0] - b[0];
  bp[1] = p[1] - b[1];
  bp[2] = p[2] - b[2];

  Scalar d3 = ab[0] * bp[0] + ab[1] * bp[1] + ab[2] * bp[2];
  Scalar d4 = ac[0] * bp[0] + ac[1] * bp[1] + ac[2] * bp[2];

  if ((d3 >= 0.0f) && (d4 <= d3)) {
    t1 = 0.0f;
    t2 = 1.0f;
    t3 = 0.0f;

    return std::sqrt((p[0] - b[0]) * (p[0] - b[0]) +
                     (p[1] - b[1]) * (p[1] - b[1]) +
                     (p[2] - b[2]) * (p[2] - b[2]));
  }

  Scalar vc = d1 * d4 - d3 * d2;

  if ((vc <= 0.0f) && (d1 >= 0.0f) && (d3 <= 0.0f)) {
    Scalar v = d1 / (d1 - d3);

    t1 = 1 - v;
    t2 = v;
    t3 = 0;

    Scalar vec[3];
    vec[0] = p[0] - (a[0] + v * ab[0]);
    vec[1] = p[1] - (a[1] + v * ab[1]);
    vec[2] = p[2] - (a[2] + v * ab[2]);

    return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  }

  Scalar cp[3];
  cp[0] = p[0] - c[0];
  cp[1] = p[1] - c[1];
  cp[2] = p[2] - c[2];

  Scalar d5 = ab[0] * cp[0] + ab[1] * cp[1] + ab[2] * cp[2];
  Scalar d6 = ac[0] * cp[0] + ac[1] * cp[1] + ac[2] * cp[2];

  if ((d6 >= 0.0f) && (d5 <= d6)) {
    t1 = 0;
    t2 = 0;
    t3 = 1;

    return std::sqrt((p[0] - c[0]) * (p[0] - c[0]) +
                     (p[1] - c[1]) * (p[1] - c[1]) +
                     (p[2] - c[2]) * (p[2] - c[2]));
  }

  Scalar vb = d5 * d2 - d1 * d6;

  if ((vb <= 0.0f) && (d2 >= 0.0f) && (d6 <= 0.0f)) {
    Scalar w = d2 / (d2 - d6);

    t1 = 1 - w;
    t2 = 0;
    t3 = w;

    Scalar vec[3];
    vec[0] = p[0] - (a[0] + w * ac[0]);
    vec[1] = p[1] - (a[1] + w * ac[1]);
    vec[2] = p[2] - (a[2] + w * ac[2]);

    return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  }

  Scalar va = d3 * d6 - d5 * d4;

  if ((va <= 0.0f) && ((d4 - d3) >= 0.0f) && ((d5 - d6) >= 0.0f)) {
    Scalar w = (d4 - d3) / ((d4 - d3) + (d5 - d6));

    t1 = 0;
    t2 = 1 - w;
    t3 = w;

    Scalar vec[3];
    vec[0] = p[0] - (b[0] + w * (c[0] - b[0]));
    vec[1] = p[1] - (b[1] + w * (c[1] - b[1]));
    vec[2] = p[2] - (b[2] + w * (c[2] - b[2]));

    return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  }

  Scalar denom = 1.0f / (va + vb + vc);
  Scalar v = vb * denom;
  Scalar w = vc * denom;
  Scalar u = 1.0 - v - w;

  t1 = u;
  t2 = v;
  t3 = w;

  Scalar vec[3];
  vec[0] = p[0] - (u * a[0] + v * b[0] + w * c[0]);
  vec[1] = p[1] - (u * a[1] + v * b[1] + w * c[1]);
  vec[2] = p[2] - (u * a[2] + v * b[2] + w * c[2]);

  return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

void print_histogram_analysis(const std::vector<VecXx> &vec, int num_bins,
                              const std::string &name, bool nonzero) {
  if (!vec.size()) return;

  VecXx bucket_mins(vec.size());
  bucket_mins.setConstant(std::numeric_limits<Scalar>::infinity());

  VecXx bucket_maxs(vec.size());
  bucket_maxs.setConstant(-std::numeric_limits<Scalar>::infinity());

  const int num_buckets = (int)vec.size();

  for_each(0, num_buckets, [&](int bucket_idx) {
    if (vec[bucket_idx].size() == 0) return;

    bucket_maxs[bucket_idx] = vec[bucket_idx].maxCoeff();
    bucket_mins[bucket_idx] = vec[bucket_idx].minCoeff();
  });

  Scalar total_min = bucket_mins.minCoeff();
  Scalar total_max = bucket_maxs.maxCoeff();
  Scalar range = (total_max - total_min);

  if (range == 0.0) return;

  std::vector<std::vector<int> > histogram_by_buckets(vec.size());

  for_each(0, num_buckets, [&](int bucket_idx) {
    const VecXx &nodes = vec[bucket_idx];
    std::vector<int> &hist = histogram_by_buckets[bucket_idx];

    hist.resize(num_bins, 0);

    const int num_nodes = (int)nodes.size();
    for (int i = 0; i < num_nodes; ++i) {
      if (nonzero && nodes[i] == 0.0) continue;

      int idx_bin = std::min(num_bins - 1, (int)((nodes[i] - total_min) /
                                                 range * (Scalar)num_bins));
      hist[idx_bin]++;
    }
  });

  std::vector<int> total_bins(num_bins, 0);
  for_each(0, num_bins, [&](int bin_idx) {
    for (int j = 0; j < num_buckets; ++j) {
      total_bins[bin_idx] += histogram_by_buckets[j][bin_idx];
    }
  });

  int total_samples = 0;
  for (int i = 0; i < num_bins; ++i) total_samples += total_bins[i];

  std::cout << "[Hist. Anal. for " << name << ", min: " << total_min
            << ", max: " << total_max << "]" << std::endl;
  for (int i = 0; i < num_bins; ++i) {
    std::cout << std::setprecision(2)
              << ((Scalar)total_bins[i] / (Scalar)total_samples * 100.0)
              << "% ";
  }
  std::cout << std::endl;
}
}  // namespace strandsim
