/**
 * \copyright 2012 Gilles Daviet
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "LevelSet.hh"

#include "../Core/ElasticStrandUtils.hh"
#include "../Render/OpenGLHeaders.hh"
#include "../Utils/MathUtilities.hh"
#include "../Utils/TextLog.hh"
#include "Bridson/levelset_util.hh"

using namespace std;

namespace strandsim {

#define MAX_LEVEL_SET_CELLS (1e6)
#define MIN_UNIQUE_INTERSECTIONS_RATIO (0.95)

template <int deg>
InterpolatedLevelSet<deg>::InterpolatedLevelSet()
    : m_dx(1.),
      m_scale(1.0),
      m_initialized(false),
      m_needsRebuilding(false),
      m_adaptativeGridSize(false),
      m_transformChanged(false),
      m_x(NULL),
      m_v(NULL),
      m_triangles(NULL),
      m_triIndices(NULL) {
  m_realRequestPositions.clear();
  m_transformedRequestPositions.clear();
  m_grad.clear();
  m_gradPosition.clear();

  m_transformMatrixAtCreation.setIdentity();
  m_currentTransformMatrix.setIdentity();
  setTransformationMatrix(m_transformMatrixAtCreation);
  setScale(1.0);
}

template <int deg>
InterpolatedLevelSet<deg>::~InterpolatedLevelSet() {}

template <int deg>
void InterpolatedLevelSet<deg>::calculateLevelSetSize(
    const FaceArray &triVertices, const Indices &triangles, const Vec3xArray &x,
    const Vec3xArray &v, Vec3x &origin, Vec3x &center) {
  Vec3x xMin(Vec3x::Zero()), xMax(Vec3x::Zero());

  if (triangles.size()) {
    // Find the mesh's bounding box
    for (int i = 0; i < 3; ++i) {
      xMin[i] = std::numeric_limits<Scalar>::infinity();
      xMax[i] = -std::numeric_limits<Scalar>::infinity();
    }

    for (unsigned t = 0; t < triangles.size(); ++t) {
      for (unsigned v = 0; v < 3; ++v) {
        const Vec3x &vertex = x[triVertices[triangles[t]](v)];

        for (unsigned int k = 0; k < 3; ++k) {
          if (vertex[k] < xMin[k]) {
            xMin[k] = vertex[k];
          }
          if (vertex[k] > xMax[k]) {
            xMax[k] = vertex[k];
          }
        }
      }
    }
  }

  center = (xMin + xMax) / 2.0;
  origin = xMin;

  //    std::cerr << "LevelSet::calculateLevelSetSize() - size:\n";
  //    std::cerr << "   origin = ( " << origin[0] << ", " << origin[1] << ", "
  //    << origin[2] << " )" << std::endl; std::cerr << "   center = ( " <<
  //    center[0] << ", " << center[1] << ", " << center[2] << " )" <<
  //    std::endl;
}

template <int deg>
void InterpolatedLevelSet<deg>::buildLevelSet(
    const FaceArray &triangles, const Indices &triIndices, const Vec3xArray &x,
    const Vec3xArray &v, const Vec3x &origin, const Scalar dx, const int nx,
    const int ny, const int nz, const Mat4x &transformMatrix) {
  // DebugStream( g_log, "" ) << "InterpolatedLevelSet<" << deg <<
  // ">::buildLevelSet";

  {
    LockGuard lock(m_levelSetMutex);
    LockGuard tlock(m_transformMutex);

    m_newOrigin = origin;
    m_newCenter = origin + Vec3x(nx, ny, nz) * dx * 0.5;

    m_x = &x;
    m_v = &v;

    // Chose a consensual naming convention
    m_triangles = &triIndices;
    m_triIndices = &triangles;

    m_dx = dx;

    m_transformMatrixAtCreation = transformMatrix;

    m_needsRebuilding = true;
    m_transformChanged = true;
    m_initialized = true;
  }

  update();  // For compatibility reasons, perform synchronous update when
             // called from figaro
}

template <int deg>
void InterpolatedLevelSet<deg>::buildLevelSetFromFile(
    std::ifstream &levelSetFile, const FaceArray &triangles,
    const Indices &triIndices, const Vec3xArray &x, const Vec3xArray &v,
    const Vec3x &origin, const Scalar dx, const int nx, const int ny,
    const int nz, const Mat4x &transformMatrix) {
  // DebugStream( g_log, "" ) << "InterpolatedLevelSet<" << deg <<
  // ">::buildLevelSet";

  {
    LockGuard lock(m_levelSetMutex);
    LockGuard tlock(m_transformMutex);

    m_newOrigin = origin;
    m_newCenter = origin + Vec3x(nx, ny, nz) * dx * 0.5;

    m_x = &x;
    m_v = &v;

    // Chose a consensual naming convention
    m_triangles = &triIndices;
    m_triIndices = &triangles;

    m_dx = dx;

    m_transformMatrixAtCreation = transformMatrix;

    m_needsRebuilding = true;
    m_transformChanged = true;
    m_initialized = true;
  }

  updateFromFile(levelSetFile);  // For compatibility reasons, perform
                                 // synchronous update when called from figaro
}

template <int deg>
void InterpolatedLevelSet<deg>::buildAdaptativeLevelSet(
    const FaceArray &triangles, const Indices &triIndices, const Vec3xArray &x,
    const Vec3xArray &v, const Mat4x &transformMatrix) {
  LockGuard lock(m_levelSetMutex);
  LockGuard tlock(m_transformMutex);

  calculateLevelSetSize(triangles, triIndices, x, v, m_newOrigin, m_newCenter);

  m_x = &x;
  m_v = &v;

  // Chose a consensual naming convention
  m_triangles = &triIndices;
  m_triIndices = &triangles;

  m_adaptativeGridSize = true;

  m_transformMatrixAtCreation = transformMatrix;

  m_needsRebuilding = true;
  m_transformChanged = true;
  m_initialized = true;
}

template <int deg>
void InterpolatedLevelSet<deg>::update() {
  {
    LockGuard lock(m_levelSetMutex);

    if (m_needsRebuilding) {
      m_needsRebuilding = false;

      m_phi.clear();
      m_phiVel.clear();

      m_origin = m_newOrigin;
      m_unpaddedOrigin = m_origin;

      if (m_origin == m_newCenter) {
        // DebugStream( g_log, "" ) << " Empty level set -- postponing
        // initialization ";
      } else {
        // DebugStream( g_log, "" ) << "Building level set with origin = " <<
        // m_origin
        //                    << " and center = " << m_newCenter;
        // DebugStream( g_log, "" ) << "Transformation matrix at creation: "
        //                    << m_transformMatrixAtCreation;

        computeLevelSet(*m_triIndices, *m_triangles, *m_x, *m_v, m_origin,
                        m_newCenter, m_dx, m_phi, m_phiVel, m_closest_tri,
                        m_adaptativeGridSize);
      }

      m_length = 2. * (m_newCenter - m_origin);

      computeAABB();
    }
  }

  {
    LockGuard lock(m_transformMutex);

    if (m_transformChanged) {
      m_transformChanged = false;
      setScale(
          1.0);  // would eventually get overridden if we did per-strand scaling
    }
  }
}

template <int deg>
void InterpolatedLevelSet<deg>::updateFromFile(std::ifstream &levelSetFile) {
  {
    LockGuard lock(m_levelSetMutex);

    if (m_needsRebuilding) {
      m_needsRebuilding = false;

      m_phi.clear();
      m_phiVel.clear();

      m_origin = m_newOrigin;
      m_unpaddedOrigin = m_origin;

      if (m_origin == m_newCenter) {
        // DebugStream( g_log, "" ) << " Empty level set -- postponing
        // initialization ";
      } else {
        // DebugStream( g_log, "" ) << "Building level set with origin = " <<
        // m_origin
        //                    << " and center = " << m_newCenter;
        // DebugStream( g_log, "" ) << "Transformation matrix at creation: "
        //                    << m_transformMatrixAtCreation;

        loadFile(levelSetFile);
      }

      m_length = 2. * (m_newCenter - m_origin);

      computeAABB();
    }
  }

  {
    LockGuard lock(m_transformMutex);

    if (m_transformChanged) {
      m_transformChanged = false;
      setScale(
          1.0);  // would eventually get overridden if we did per-strand scaling
    }
  }
}

template <int deg>
void InterpolatedLevelSet<deg>::setTransformationMatrix(const Mat4x &i_matrix) {
  LockGuard lock(m_transformMutex);

  m_currentTransformMatrix = i_matrix;

  m_transformChanged = true;
}

template <int deg>
void InterpolatedLevelSet<deg>::setScale(const Scalar scale) {
  m_scale = scale;
  Mat4x scaling;
  const Vec3x &center = m_origin + m_length * .5;
  scaling << m_scale, 0., 0., (1. - m_scale) * center[0], 0., m_scale, 0.,
      (1. - m_scale) * center[1], 0., 0., m_scale, (1. - m_scale) * center[2],
      0., 0., 0., 1.;

  m_relativeTransformMatrix = m_transformMatrixAtCreation *
                              (scaling * m_currentTransformMatrix).inverse();
  m_invRelativeTransformMatrix = m_relativeTransformMatrix.inverse();

  computeAABB();
}

// For debug purposes only
template <int deg>
void InterpolatedLevelSet<deg>::set(bridson::Array3x phi, Vec3x origin,
                                    Vec3x center, Scalar dx) {
  m_phi = phi;
  m_phiVel.resize(phi.ni, phi.nj, phi.nk);
  m_origin = origin;
  m_length = 2. * (center - origin);
  m_dx = dx;
  m_scale = 1.f;
  m_relativeTransformMatrix.setIdentity();
  m_invRelativeTransformMatrix.setIdentity();
  m_initialized = true;
}

template <int deg>
void InterpolatedLevelSet<deg>::computeAABB() {
  Vec3x worldOrig = transformPoint(m_invRelativeTransformMatrix, m_origin);

  Vec3x locEdge;
  Mat3x worldEdges;

  for (unsigned k = 0; k < 3; ++k) {
    locEdge.setZero();
    locEdge(k) = m_length(k);
    worldEdges.col(k) =
        (m_invRelativeTransformMatrix.block<3, 3>(0, 0) * locEdge);
  }

  m_aaBBMin = worldOrig;
  m_aaBBMax = worldOrig;

  Vec3x vertex;
  for (unsigned i = 1; i < 8; ++i) {
    vertex = worldOrig;
    for (unsigned k = 0; k < 3; ++k) {
      vertex += ((i & (1 << k)) >> k) * worldEdges.col(k);
    }

    m_aaBBMin = m_aaBBMin.cwiseMin(vertex);
    m_aaBBMax = m_aaBBMax.cwiseMax(vertex);
  }

  //    std::cout << " AABB min : " << m_aaBBMin.transpose() << std::endl ;
  //    std::cout << " AABB max : " << m_aaBBMax.transpose() << std::endl ;
}

template <int deg>
bool InterpolatedLevelSet<deg>::intersectsAABB(
    const std::pair<Vec3x, Vec3x> &aaBB) const {
  const Vec3x &min = aaBB.first;
  const Vec3x &max = aaBB.second;

  for (unsigned k = 0; k < 3; ++k) {
    if (max[k] < m_aaBBMin[k] || m_aaBBMax[k] < min[k]) return false;
  }

  return true;
}

template <int deg>
Vec3x InterpolatedLevelSet<deg>::getWorldSpacePos(int i, int j, int k) const {
  Vec4x worPos = m_invRelativeTransformMatrix *
                 Vec4x(i * m_dx + m_origin[0], j * m_dx + m_origin[1],
                       k * m_dx + m_origin[2], 1.0f);
  return Vec3x(worPos[0], worPos[1], worPos[2]);
}

template <int deg>
Scalar InterpolatedLevelSet<deg>::getLevelSetValue(int i, int j, int k) const {
  return getInterpolated(i, j, k, 0, 0, 0) * m_scale;
}

template <int deg>
Scalar InterpolatedLevelSet<deg>::getLevelSetValue(const Vec3x &x) const {
  assert(m_initialized);

  // Compute the evaluation point in the level set's original coordinates
  const Vec4x &samplePoint =
      m_relativeTransformMatrix * Vec4x(x[0], x[1], x[2], 1.0f);
  // Compute the grid coordinates, integer and floating part
  int i, j, k;

  Scalar fi = (samplePoint[0] - m_origin[0]) / m_dx;
  Scalar fj = (samplePoint[1] - m_origin[1]) / m_dx;
  Scalar fk = (samplePoint[2] - m_origin[2]) / m_dx;

  static const int bound = (deg - 1) / 2;

  bridson::get_barycentric(fi, i, fi, bound, m_phi.ni - bound);
  bridson::get_barycentric(fj, j, fj, bound, m_phi.nj - bound);
  bridson::get_barycentric(fk, k, fk, bound, m_phi.nk - bound);

  Scalar dist = getInterpolated(i, j, k, fi, fj, fk);

  // Scale the distance back into world coordinates
  return m_scale * dist;
}

template <int deg>
void InterpolatedLevelSet<deg>::getGradient(const Vec3x &x, Vec3x &grad,
                                            bool rotated) const {
  assert(m_initialized);

  // Compute the evaluation point in the level set's original coordinates
  const Vec4x &samplePoint =
      m_relativeTransformMatrix * Vec4x(x[0], x[1], x[2], 1.0f);

  // Compute the grid coordinates, integer and floating part
  int i, j, k;

  Scalar fi = (samplePoint[0] - m_origin[0]) / m_dx;
  Scalar fj = (samplePoint[1] - m_origin[1]) / m_dx;
  Scalar fk = (samplePoint[2] - m_origin[2]) / m_dx;

  bridson::get_barycentric(fi, i, fi, 1, m_phi.ni - 1);
  bridson::get_barycentric(fj, j, fj, 1, m_phi.nj - 1);
  bridson::get_barycentric(fk, k, fk, 1, m_phi.nk - 1);

  // If needed, transform the gradient vector back into world coordinates, but
  // without translation
  if (rotated) {
    Mat3x M = m_invRelativeTransformMatrix.block<3, 3>(0, 0);
    M = m_scale * M;
    grad = M * getInterpolatedGradient(i, j, k, fi, fj, fk);
  } else {
    grad = getInterpolatedGradient(i, j, k, fi, fj, fk);
  }
}

template <int deg>
void InterpolatedLevelSet<deg>::getHessian(const Vec3x &x, Mat3x &hess,
                                           bool rotated) const {
  assert(m_initialized);

  // Compute the evaluation point in the level set's original coordinates
  const Vec4x &samplePoint =
      m_relativeTransformMatrix * Vec4x(x[0], x[1], x[2], 1.0f);

  // Compute the grid coordinates, integer and floating part
  int i, j, k;

  Scalar fi = (samplePoint[0] - m_origin[0]) / m_dx;
  Scalar fj = (samplePoint[1] - m_origin[1]) / m_dx;
  Scalar fk = (samplePoint[2] - m_origin[2]) / m_dx;

  bridson::get_barycentric(fi, i, fi, 1, m_phi.ni - 1);
  bridson::get_barycentric(fj, j, fj, 1, m_phi.nj - 1);
  bridson::get_barycentric(fk, k, fk, 1, m_phi.nk - 1);

  hess = getInterpolatedHessian(i, j, k, fi, fj, fk);

  if (rotated) {
    hess = m_scale * m_invRelativeTransformMatrix.block<3, 3>(0, 0) * hess *
           m_relativeTransformMatrix.block<3, 3>(0, 0);
  }
}

template <int deg>
Scalar InterpolatedLevelSet<deg>::getEdgeValue(const Vec3x &start,
                                               const Vec3x &end,
                                               Coord &min) const {
  assert(m_initialized);

  // Compute the evaluation point in the level set's original coordinates
  Vec4x v0 =
      (m_relativeTransformMatrix * Vec4x(start[0], start[1], start[2], 1.0f));
  Vec4x v1 = (m_relativeTransformMatrix * Vec4x(end[0], end[1], end[2], 1.0f));

  const Vec3x &sPoint = v0.segment<3>(0);
  const Vec3x &ePoint = v1.segment<3>(0);

  getEdgeSamplingPoints(sPoint, ePoint, min);

  // Scale the distance back into world coordinates
  return m_scale * min.d;
}

template <int deg>
void InterpolatedLevelSet<deg>::getInteriorPoints(const Vec3x &start,
                                                  const Vec3x &end,
                                                  const Scalar thickness,
                                                  std::vector<Coord> &points,
                                                  bool stopAtFirstMin) const {
  assert(m_initialized);

  // Compute the evaluation point in the level set's original coordinates
  Vec4x v0 =
      (m_relativeTransformMatrix * Vec4x(start[0], start[1], start[2], 1.0f));

  const Vec3x &sPoint = v0.segment<3>(0);

  Vec4x v1 = (m_relativeTransformMatrix * Vec4x(end[0], end[1], end[2], 1.0f));
  const Vec3x &ePoint = v1.segment<3>(0);

  Coord min;
  getEdgeSamplingPoints(sPoint, ePoint, min, &points, thickness / m_scale,
                        stopAtFirstMin);
}

template <int deg>
void InterpolatedLevelSet<deg>::getInteriorIntervals(
    const Vec3x &start, const Vec3x &end, const Scalar thickness,
    std::vector<Interval> &intervals) const {
  assert(m_initialized);

  // Compute the evaluation point in the level set's original coordinates
  Vec4x v0 =
      (m_relativeTransformMatrix * Vec4x(start[0], start[1], start[2], 1.0f));

  const Vec3x &sPoint = v0.segment<3>(0);

  Vec4x v1 = (m_relativeTransformMatrix * Vec4x(end[0], end[1], end[2], 1.0f));
  const Vec3x &ePoint = v1.segment<3>(0);

  getEdgeIntervals(sPoint, ePoint, thickness / m_scale, intervals);
}

template <int deg>
void InterpolatedLevelSet<deg>::getEdgeGradient(const Vec3x &start,
                                                const Vec3x &end,
                                                const Coord &min, Vec3x &grad,
                                                int *gradLambdaAxis) const {
  const Vec3x &colPoint = (1. - min.t) * start + min.t * end;

  Vec3x dLS_dP;
  getGradient(colPoint, dLS_dP, false);

  Mat3x dP_dPi = Mat3x::Identity();

  int gradAxis = -1;
  if (!isSmall(min.t) && !isSmall(1. - min.t)) {
    // Find closest plane ( which normal is gradLambda axis)
    Scalar planeDist = 1.;
    for (unsigned k = 0; k < 3; ++k) {
      const Scalar r = std::min(min.l[k], 1 - min.l[k]);
      if (r < planeDist) {
        gradAxis = k;
        planeDist = r;
      }
    }

    const Vec3x &dir =
        (m_relativeTransformMatrix.block<3, 3>(0, 0) * (end - start));

    const Scalar dGradLambda = -1. / dir[gradAxis];

    // dP_dPi += dir * gradAlpha.transpose() ;
    dP_dPi.col(gradAxis) += dir * dGradLambda;
  }

  if (gradLambdaAxis) *gradLambdaAxis = gradAxis;

  grad = m_scale * m_invRelativeTransformMatrix.block<3, 3>(0, 0) *
         dP_dPi.transpose() * dLS_dP;
}

template <int deg>
void InterpolatedLevelSet<deg>::getEdgeHessian(
    const Vec3x &start, const Vec3x &end, const Scalar lambda,
    const int gradLambdaAxis, Eigen::Matrix<Scalar, 6, 6> &hess) const {
  const bool interior = gradLambdaAxis >= 0;

  const Vec3x &colPoint = (1. - lambda) * start + lambda * end;

  Vec3x dLS_dP;
  getGradient(colPoint, dLS_dP, false);

  const Vec3x &dir =
      (m_relativeTransformMatrix.block<3, 3>(0, 0) * (end - start));

  const Scalar dGradLambda = interior ? -1. / dir[gradLambdaAxis] : 0.;
  const Scalar dHessLambda = dGradLambda * dGradLambda;

  Mat3x d2LS_dP2;
  getHessian(colPoint, d2LS_dP2, false);

  Mat3x dP_dPi = Mat3x::Identity();
  if (interior) {
    dP_dPi.col(gradLambdaAxis) += dir * dGradLambda;
  }
  const Mat3x hessLS = dP_dPi.transpose() * d2LS_dP2 * dP_dPi;

  Mat3x hessP;
  hessP.setZero();

  for (unsigned i = 0; i < 2; ++i) {
    const Scalar lambdai = i ? lambda : (lambda - 1);

    for (unsigned j = 0; j < 2; ++j) {
      const Scalar lambdaj = j ? lambda : (lambda - 1);
      const int sign = i == j ? 1 : -1;

      if (interior) {
        hessP.setZero();

        for (unsigned k = 0; k < 3; ++k) {
          // Mat3x dPk_dPi ;
          // Vec3x dPik_dPi  = Vec3x::Zero() ;

          // dPk_dPi = lambdai * dPik_dPi * gradLambda.transpose() ;
          hessP(k, gradLambdaAxis) += dLS_dP[k] * lambdai * dGradLambda;
          // dPk_dPi += lambdaj * gradLambda * dPik_dPi.transpose() ;
          hessP(gradLambdaAxis, k) += dLS_dP[k] * lambdaj * dGradLambda;
          // dPk_dPi += (lambdai + lambdaj) * dir[ k ] * hessLambda ;
          hessP(gradLambdaAxis, gradLambdaAxis) +=
              dLS_dP[k] * (lambdai + lambdaj) * dir[k] * dHessLambda;

          // hessP +=  dLS_dP[k] * dPk_dPi ;
        }
      }

      hess.block<3, 3>(3 * j, 3 * i) =
          m_invRelativeTransformMatrix.block<3, 3>(0, 0) * m_scale * sign *
          (lambdai * lambdaj * hessLS + hessP) *
          m_relativeTransformMatrix.block<3, 3>(0, 0);
    }
  }
  // std::cout << "--\n " << hess << std::endl ;
}

template <int deg>
Scalar InterpolatedLevelSet<deg>::getIntegralValue(const Interval &iv) const {
  std::cerr << "Not implemented" << std::endl;
  return -1;
}

template <int deg>
void InterpolatedLevelSet<deg>::getIntegralGradient(const Interval &iv,
                                                    Vec3x &gradStart,
                                                    Vec3x &gradEnd) const {
  std::cerr << "Not implemented" << std::endl;
}

template <int deg>
void InterpolatedLevelSet<deg>::getIntegralHessian(const Interval &iv,
                                                   Mat6x &hess) const {
  std::cerr << "Not implemented" << std::endl;
}

template <>
Scalar InterpolatedLevelSet<1>::getIntegralValue(const Interval &iv) const {
  Vec4x coeffs;
  computeTrilinearValueCoeffs(m_phi, iv, coeffs);

  const Vec4x intTi((iv.second.t - iv.first.t),
                    (std::pow(iv.second.t, 2) - std::pow(iv.first.t, 2)) / 2.,
                    (std::pow(iv.second.t, 3) - std::pow(iv.first.t, 3)) / 3.,
                    (std::pow(iv.second.t, 4) - std::pow(iv.first.t, 4)) / 4.);

  return m_scale * intTi.dot(coeffs);
}

template <>
void InterpolatedLevelSet<1>::getIntegralGradient(const Interval &iv,
                                                  Vec3x &gradStart,
                                                  Vec3x &gradEnd) const {
  Mat3x coeffs;
  computeBilinearGradientCoeffs(m_phi, iv, coeffs);

  coeffs =
      m_scale / m_dx * m_invRelativeTransformMatrix.block<3, 3>(0, 0) * coeffs;

  const Vec4x intTi((iv.second.t - iv.first.t),
                    (std::pow(iv.second.t, 2) - std::pow(iv.first.t, 2)) / 2.,
                    (std::pow(iv.second.t, 3) - std::pow(iv.first.t, 3)) / 3.,
                    (std::pow(iv.second.t, 4) - std::pow(iv.first.t, 4)) / 4.);

  gradEnd = coeffs * intTi.segment<3>(1);
  gradStart = coeffs * (intTi.segment<3>(0) - intTi.segment<3>(1));

  // Apporimate calculation ; alpha outisde integral
  //    const Scalar alpha = .5 * ( iv.second.t + iv.first.t ) ;
  //    const Vec3x grad = coeffs * ( intTi.segment< 3 > ( 0 ) ) ;
  //    gradEnd = alpha * grad ;
  //    gradStart = (1 - alpha ) * grad ;
}

template <>
void InterpolatedLevelSet<1>::getIntegralHessian(const Interval &iv,
                                                 Mat6x &hess) const {
  Mat3x HA, HB;
  computeLinearHessianCoeffs(m_phi, iv, HA, HB);

  HA = m_scale / (m_dx * m_dx) *
       m_invRelativeTransformMatrix.block<3, 3>(0, 0) * HA *
       m_relativeTransformMatrix.block<3, 3>(0, 0);
  HB = m_scale / (m_dx * m_dx) *
       m_invRelativeTransformMatrix.block<3, 3>(0, 0) * HB *
       m_relativeTransformMatrix.block<3, 3>(0, 0);

  const Vec4x intTi((iv.second.t - iv.first.t),
                    (std::pow(iv.second.t, 2) - std::pow(iv.first.t, 2)) / 2.,
                    (std::pow(iv.second.t, 3) - std::pow(iv.first.t, 3)) / 3.,
                    (std::pow(iv.second.t, 4) - std::pow(iv.first.t, 4)) / 4.);

  const Mat3x intHt = HA * intTi[2] + HB * intTi[1];

  hess.block<3, 3>(3, 3) = HA * intTi[3] + HB * intTi[2];
  hess.block<3, 3>(3, 0) = intHt - hess.block<3, 3>(3, 3);
  hess.block<3, 3>(0, 3) = hess.block<3, 3>(3, 0);
  hess.block<3, 3>(0, 0) =
      hess.block<3, 3>(3, 3) + HA * intTi[1] + HB * intTi[0] - 2 * (intHt);

  // Apporimate calculation ; alpha outisde integral
  //    const Scalar alpha = .5 * ( iv.second.t + iv.first.t ) ;
  //    const Mat3x J = ( HA * intTi[ 1 ] + HB * intTi[ 0 ] ) ;
  //    hess.block< 3, 3 > ( 3, 3 ) = alpha * alpha * J ;
  //    hess.block< 3, 3 > ( 3, 0 ) = alpha * (1 - alpha ) * J ;
  //    hess.block< 3, 3 > ( 0, 3 ) = (1 - alpha ) * alpha * J ;
  //    hess.block< 3, 3 > ( 0, 0 ) = (1 - alpha ) * (1 - alpha ) * J ;

  // std::cout << hess << std::endl ;
}

template <int deg>
void InterpolatedLevelSet<deg>::getEdgeIntervals(
    const Vec3x &start, const Vec3x &end, const Scalar maxVal,
    std::vector<Interval> &intervals) const {
  std::vector<Coord> points;
  Coord min;
  getEdgeSamplingPoints(start, end, min, &points, maxVal, false);

  if (min.d > maxVal) return;

  // FIXME for now we use a linear approximation to compute intersection with
  // "thickness" isosurface We don't need a high precision for this (do we ?),
  // so I guess it should be ok

  Interval iv;
  for (unsigned i = 1; i < points.size(); ++i) {
    iv.first = points[i - 1];
    iv.second = points[i];

    if (iv.second.d > maxVal) {
      iv.second.t =
          iv.first.t + (iv.second.t - iv.first.t) *
                           ((maxVal - iv.first.d) / (iv.second.d - iv.first.d));
      iv.second.f = (1. - iv.second.t) * iv.first.f + iv.second.t * iv.second.f;

      computeBarycentric(iv.second);
      // iv.second.d = getInterpolated( iv.second );

    } else if (iv.first.d > maxVal) {
      iv.first.t =
          iv.first.t + (iv.second.t - iv.first.t) *
                           ((maxVal - iv.first.d) / (iv.second.d - iv.first.d));
      iv.first.f = (1. - iv.first.t) * iv.first.f + iv.first.t * iv.second.f;

      computeBarycentric(iv.first);
      // iv.first.d = getInterpolated( iv.first );
    }
    intervals.push_back(iv);
    normalizeInterval(intervals.back());
  }
}

template <int deg>
Scalar InterpolatedLevelSet<deg>::getLevelSetValueVelocity(const Vec3x &x,
                                                           Vec3x &v) const {
  Vec4x samplePoint(x[0], x[1], x[2], 1);

  samplePoint = m_relativeTransformMatrix * samplePoint;

  int i, j, k;

  Scalar fi = ((Scalar)samplePoint[0] - m_origin[0]) / m_dx;
  Scalar fj = ((Scalar)samplePoint[1] - m_origin[1]) / m_dx;
  Scalar fk = ((Scalar)samplePoint[2] - m_origin[2]) / m_dx;

  static const int bound = (deg - 1) / 2;

  bridson::get_barycentric(fi, i, fi, bound, m_phi.ni - bound);
  bridson::get_barycentric(fj, j, fj, bound, m_phi.nj - bound);
  bridson::get_barycentric(fk, k, fk, bound, m_phi.nk - bound);

  Scalar dist = bridson::trilerp(
      m_phi(i, j, k), m_phi(i + 1, j, k), m_phi(i, j + 1, k),
      m_phi(i + 1, j + 1, k), m_phi(i, j, k + 1), m_phi(i + 1, j, k + 1),
      m_phi(i, j + 1, k + 1), m_phi(i + 1, j + 1, k + 1), fi, fj, fk);

  v = bridson::trilerp(m_phiVel(i, j, k), m_phiVel(i + 1, j, k),
                       m_phiVel(i, j + 1, k), m_phiVel(i + 1, j + 1, k),
                       m_phiVel(i, j, k + 1), m_phiVel(i + 1, j, k + 1),
                       m_phiVel(i, j + 1, k + 1), m_phiVel(i + 1, j + 1, k + 1),
                       fi, fj, fk);

  return m_scale * dist;
}

template <int deg>
void InterpolatedLevelSet<deg>::getEdgeSamplingPoints(
    const Vec3x &start, const Vec3x &end, Coord &min,
    std::vector<Coord> *points, Scalar maxVal, bool stopAtFirstMin) const {
  static const Scalar eps = SMALL_NUMBER<Scalar>();

  min.d = std::numeric_limits<Scalar>::infinity();

  const Scalar invCellSize = 1.f / m_dx;

  const Vec3x &fs = (start - m_origin) * invCellSize;
  const Vec3x &dP = (end - start) * invCellSize;
  Vec3x invdP;

  Vec3x fwd;  // fwd[i] = 1 if dP[i].e_i > 0

  // First find the first and last t that are inside the level set grid
  bool notColliding = false;
  Scalar tin = 0.f, tout = 1.f;
  for (unsigned k = 0; k < 3; ++k) {
    bool rev = dP[k] < 0;
    fwd[k] = !rev;

    if (isSmall(dP[k])) {
      invdP[k] = (2 * fwd[k] - 1) / (eps);  // relatively close to infinity
    } else {
      invdP[k] = 1. / dP[k];
    }

    // Account for the projection in compute_barycentric for points < 1 or > nx
    // - 2
    Scalar P0 = (1 - fs[k]) * invdP[k];
    Scalar P1 = P0 + (m_length[k] * invCellSize - 3) * invdP[k];

    if (rev) std::swap(P0, P1);

    if (P0 > 1 || P1 < 0) {
      notColliding = true;
      break;
    }

    tin = std::max(tin, P0);
    tout = std::min(tout, P1);
  }

  if (notColliding || tin > tout || invdP.lpNorm<Eigen::Infinity>() < eps)
    return;

  // Then find all intersections with grid
  Coord cur, prev;
  cur.t = tin;

  Vec3x dx;

  bool done = false;
  bool inside = false;
  bool pushed_back = true;

  bool decreasing = false;

  while (true) {
    if (tout - cur.t < eps) {
      cur.t = tout;
      done = true;
    }

    cur.f = fs + cur.t * dP;

    computeBarycentric(cur);
    cur.d = getTrilerp(cur);

    if (cur.d < min.d) {
      decreasing = true;
      min = cur;
    } else if (stopAtFirstMin && decreasing) {
      // break ;
    }

    if (points) {
      if (stopAtFirstMin) {
        if (decreasing) (*points).push_back(cur);
      } else {
        // Also add before and after points ( required by getInteriorIntervals )
        const bool wasInside = inside;
        inside = cur.d <= maxVal;

        if (wasInside || inside) {
          if (!pushed_back) {
            (*points).push_back(prev);
          }

          (*points).push_back(cur);
          pushed_back = true;

        } else {
          pushed_back = false;
          prev = cur;
        }
      }
    }

    if (done) break;

    // Increment to next intersecting plane
    dx = (fwd - cur.l).cwiseProduct(invdP);

    // Avoid getting stuck by numerical imprecisions
    Scalar dt = 1.;
    for (int k = 0; k < 3; ++k) {
      if (dx[k] < eps) {
        dx[k] += std::fabs(invdP[k]);  // Jump to the next cell
      }
      if (dx[k] < dt) dt = dx[k];
    }

    cur.t += dt;
  }
}

template <int deg>
unsigned InterpolatedLevelSet<deg>::getClosestTriangle(const Vec3x &x,
                                                       Scalar &t1, Scalar &t2,
                                                       Scalar &t3) const {
  assert(m_initialized);

  // Compute the evaluation point in the level set's original coordinates
  const Vec4x &samplePoint =
      m_relativeTransformMatrix * Vec4x(x[0], x[1], x[2], 1.0f);

  Scalar fi = ((Scalar)samplePoint[0] - m_origin[0]) / m_dx;
  Scalar fj = ((Scalar)samplePoint[1] - m_origin[1]) / m_dx;
  Scalar fk = ((Scalar)samplePoint[2] - m_origin[2]) / m_dx;

  int i, j, k;
  bridson::get_barycentric(fi, i, fi, 0, m_phi.ni);
  bridson::get_barycentric(fj, j, fj, 0, m_phi.nj);
  bridson::get_barycentric(fk, k, fk, 0, m_phi.nk);

  const unsigned t = m_closest_tri(i, j, k);

  if (t == -1U) return -1U;

  const unsigned a = (*m_triIndices)[t](0);
  const unsigned b = (*m_triIndices)[t](1);
  const unsigned c = (*m_triIndices)[t](2);

  const Vec3x P(samplePoint.segment<3>(0));

  const Vec3x &A = (*m_x)[a];
  const Vec3x &B = (*m_x)[b];
  const Vec3x &C = (*m_x)[c];

  // Assume triangle is not degenerate ; if it is, we should not be here, right
  // ?

  const Vec3x &e1 = (B - A).normalized();
  const Vec3x &e2 = (C - A - e1 * (C - A).dot(e1)).normalized();
  const Vec3x &n = (e1.cross(e2)).normalized();

  t1 = (P - A).dot(e1);
  t2 = (P - A).dot(e2);
  t3 = (P - A).dot(n);

  return t;
}

template <int deg>
Vec3x InterpolatedLevelSet<deg>::getCurrentPosition(unsigned t, Scalar t1,
                                                    Scalar t2,
                                                    Scalar t3) const {
  const unsigned a = (*m_triIndices)[t](0);
  const unsigned b = (*m_triIndices)[t](1);
  const unsigned c = (*m_triIndices)[t](2);

  const Vec3x &A = (*m_x)[a];
  const Vec3x &B = (*m_x)[b];
  const Vec3x &C = (*m_x)[c];

  const Vec3x &e1 = (B - A).normalized();
  const Vec3x &e2 = (C - A - e1 * (C - A).dot(e1)).normalized();
  const Vec3x &n = (e1.cross(e2)).normalized();

  const Vec3x P(A + e1 * t1 + e2 * t2 + n * t3);

  const Vec4x eiP(P[0], P[1], P[2], 1.f);

  const Vec4x &worlPos = m_invRelativeTransformMatrix * eiP;

  return worlPos.segment<3>(0);
}

template <int deg>
void InterpolatedLevelSet<deg>::drawInterior(const Scalar scale, Scalar center,
                                             Scalar thickness) {
  if (m_phi.empty()) return;

  unsigned vId = 0;
  std::vector<GLfloat> points;
  std::vector<GLfloat> colors;

  {
    LockGuard lock(m_levelSetMutex);

    const std::vector<GLfloat>::size_type num_points =
        3 * (m_phi.ni - 1) * (m_phi.nj - 1) * (m_phi.nk - 1);
    assert(num_points < points.max_size());
    points.resize(num_points);
    assert(num_points < colors.max_size());
    colors.resize(num_points);

    Mat4x transform;
    {
      LockGuard lock(m_transformMutex);
      transform = m_invRelativeTransformMatrix;
    }

    for (int i = 0; i < m_phi.ni - 1; ++i) {
      for (int j = 0; j < m_phi.nj - 1; ++j) {
        for (int k = 0; k < m_phi.nk - 1; ++k) {
          if (m_phi(i, j, k) <= (center + thickness) &&
              m_phi(i, j, k) >= (center - thickness)) {
            const Vec3f &color = Vec3f(1.0f, 0.f, 0.f);

            Vec3f::Map(&points[3 * vId]) =
                (transform * Vec4x(m_origin[0] + i * m_dx,
                                   m_origin[1] + j * m_dx,
                                   m_origin[2] + k * m_dx, 1.0))
                    .template segment<3>(0)
                    .cast<float>();
            Vec3f::Map(&colors[3 * vId]) = color;

            ++vId;
          }
        }
      }
    }
  }

  glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
  glPushAttrib(GL_COLOR);

  glDisable(GL_LIGHTING);
  glPointSize(3.0);

  glEnableClientState(GL_COLOR_ARRAY);
  glColorPointer(3, GL_FLOAT, 0, &colors[0]);
  // glColor3f( 0.f, 0.f, 1.f) ;

  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, &points[0]);

  glDrawArrays(GL_POINTS, 0, vId);

  // deactivate vertex arrays after drawing
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);

  glPopAttrib();
  glPopAttrib();
}

template <int deg>
void InterpolatedLevelSet<deg>::draw(const Scalar scale) const {
  if (m_phi.empty()) return;

  unsigned vId = 0;
  std::vector<GLfloat> points;
  std::vector<GLfloat> colors;

  {
    LockGuard lock(m_levelSetMutex);

    const std::vector<GLfloat>::size_type num_points =
        6 * (m_phi.ni - 1) * (m_phi.nj - 1) * (m_phi.nk - 1);
    assert(num_points < points.max_size());
    points.resize(num_points);
    assert(num_points < colors.max_size());
    colors.resize(num_points);

    Mat4x transform;
    {
      LockGuard lock(m_transformMutex);
      transform = m_invRelativeTransformMatrix;
    }

    Vec3x gradient;
    for (int i = 0; i < m_phi.ni - 1; ++i) {
      for (int j = 0; j < m_phi.nj - 1; ++j) {
        for (int k = 0; k < m_phi.nk - 1; ++k) {
          Vec3f color;
          if (m_phi(i, j, k) < 0) {
            color = Vec3f(0.f, 0.f, 1.f);
          } else {
            color = Vec3f(1.0, m_phi(i, j, k) / 10.0, 1.0);
          }

          const Vec3x &pointA =
              (transform * Vec4x(m_origin[0] + i * m_dx, m_origin[1] + j * m_dx,
                                 m_origin[2] + k * m_dx, 1.0))
                  .template segment<3>(0);

          Vec3f::Map(&points[3 * vId]) = pointA.cast<float>();
          Vec3f::Map(&colors[3 * vId]) = .5 * color;

          ++vId;

          getGradient(pointA, gradient);

          Vec3f::Map(&points[3 * vId]) =
              (pointA + scale * gradient).cast<float>();
          Vec3f::Map(&colors[3 * vId]) = color;

          ++vId;
        }
      }
    }
  }

  glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT);
  glPushAttrib(GL_COLOR);

  glDisable(GL_LIGHTING);
  glLineWidth(1);

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);

  glVertexPointer(3, GL_FLOAT, 0, &points[0]);
  glColorPointer(3, GL_FLOAT, 0, &colors[0]);

  glDrawArrays(GL_LINES, 0, points.size() / 3);

  // deactivate vertex arrays after drawing
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);

  glPopAttrib();
  glPopAttrib();
}

template <int deg>
void InterpolatedLevelSet<deg>::writeFileAsync(const std::string &szfn) const {
  AsyncDataPack *data = new AsyncDataPack();
  data->m_szfn = szfn;
  data->m_origin = m_origin;
  data->m_dx = m_dx;
  data->m_phi = m_phi;
  data->m_phiVel = m_phiVel;

  std::thread t(std::bind(
      [](AsyncDataPack *data) {
        std::ofstream levelSetFile(data->m_szfn.c_str(), std::ios::binary);

        levelSetFile.write(reinterpret_cast<const char *>(&(data->m_origin[0])),
                           sizeof(Scalar) * 3);
        levelSetFile.write(reinterpret_cast<const char *>(&(data->m_dx)),
                           sizeof(Scalar));
        levelSetFile.write(reinterpret_cast<const char *>(&(data->m_phi.ni)),
                           sizeof(int));
        levelSetFile.write(reinterpret_cast<const char *>(&(data->m_phi.nj)),
                           sizeof(int));
        levelSetFile.write(reinterpret_cast<const char *>(&(data->m_phi.nk)),
                           sizeof(int));
        levelSetFile.write(
            reinterpret_cast<const char *>(&(data->m_phi(0, 0, 0))),
            sizeof(Scalar) * data->m_phi.ni * data->m_phi.nj * data->m_phi.nk);

        for (int i = 0; i < data->m_phi.ni; ++i)
          for (int j = 0; j < data->m_phi.nj; ++j)
            for (int k = 0; k < data->m_phi.nk; ++k)
              levelSetFile.write(
                  reinterpret_cast<const char *>(&(data->m_phiVel(i, j, k)[0])),
                  sizeof(Scalar) * 3);

        levelSetFile.flush();
        levelSetFile.close();

        std::cout << "[written to cache " << data->m_szfn << "]" << std::endl;
        delete data;
      },
      data));

  t.detach();
}

template <int deg>
void InterpolatedLevelSet<deg>::writeFile(std::ofstream &levelSetFile) const {
  levelSetFile.write(reinterpret_cast<const char *>(&m_origin[0]),
                     sizeof(Scalar) * 3);
  levelSetFile.write(reinterpret_cast<const char *>(&m_dx), sizeof(Scalar));
  levelSetFile.write(reinterpret_cast<const char *>(&m_phi.ni), sizeof(int));
  levelSetFile.write(reinterpret_cast<const char *>(&m_phi.nj), sizeof(int));
  levelSetFile.write(reinterpret_cast<const char *>(&m_phi.nk), sizeof(int));
  levelSetFile.write(reinterpret_cast<const char *>(&m_phi(0, 0, 0)),
                     sizeof(Scalar) * m_phi.ni * m_phi.nj * m_phi.nk);

  for (int i = 0; i < m_phi.ni; ++i)
    for (int j = 0; j < m_phi.nj; ++j)
      for (int k = 0; k < m_phi.nk; ++k)
        levelSetFile.write(
            reinterpret_cast<const char *>(&m_phiVel(i, j, k)[0]),
            sizeof(Scalar) * 3);
}

template <int deg>
void InterpolatedLevelSet<deg>::loadFile(std::ifstream &levelSetFile) {
  levelSetFile.read(reinterpret_cast<char *>(&m_origin[0]), sizeof(Scalar) * 3);
  levelSetFile.read(reinterpret_cast<char *>(&m_dx), sizeof(Scalar));

  int nx, ny, nz;
  levelSetFile.read(reinterpret_cast<char *>(&nx), sizeof(int));
  levelSetFile.read(reinterpret_cast<char *>(&ny), sizeof(int));
  levelSetFile.read(reinterpret_cast<char *>(&nz), sizeof(int));

  m_phi.resize(nx, ny, nz);
  levelSetFile.read(reinterpret_cast<char *>(&m_phi(0, 0, 0)),
                    sizeof(Scalar) * m_phi.ni * m_phi.nj * m_phi.nk);

  m_phiVel.resize(nx, ny, nz);
  for (int i = 0; i < m_phi.ni; ++i)
    for (int j = 0; j < m_phi.nj; ++j)
      for (int k = 0; k < m_phi.nk; ++k)
        levelSetFile.read(reinterpret_cast<char *>(&m_phiVel(i, j, k)[0]),
                          sizeof(Scalar) * 3);
}

template <int deg>
size_t InterpolatedLevelSet<deg>::getFileSize(const Vec3x &origin,
                                              const Vec3x &center,
                                              const Scalar &dx) const {
  const int exact_band = 1;

  const Vec3x exactOrigin = origin;  // origin will be padded later
  const Vec3x length = 2 * (center - exactOrigin);
  const Scalar minDxCube =
      length[0] * length[1] * length[2] / MAX_LEVEL_SET_CELLS;

  int ni, nj, nk;

  // Padding (to get nice zeros on border)
  ni = std::ceil(length[0] / dx) + 6;
  nj = std::ceil(length[1] / dx) + 6;
  nk = std::ceil(length[2] / dx) + 6;

  return sizeof(Scalar) * 3 + sizeof(Scalar) + sizeof(int) * 3 +
         sizeof(Scalar) * ni * nj * nk * 4;
}

template <int deg>
void InterpolatedLevelSet<deg>::computeBarycentric(Coord &c) const {
  const Eigen::Vector3i clow(1, 1, 1);
  const Eigen::Vector3i chigh(m_phi.ni - 1, m_phi.nj - 1, m_phi.nk - 1);

  bridson::get_barycentric(c.f, c.c, c.l, clow, chigh);
}

template <int deg>
void InterpolatedLevelSet<deg>::normalizeInterval(Interval &i) const {
  // Ensures that the two barycentric coordinates of an interval point to the
  // same cell
  for (unsigned k = 0; k < 3; ++k) {
    if (i.first.c[k] < i.second.c[k]) {
      if (i.first.c[k] + 2 == i.second.c[k]) {
        i.first.l[k] = 0.;
        i.first.c[k]++;
        i.second.l[k] = 1.;
        i.second.c[k]--;
      } else if ((1. - i.first.l[k]) < i.second.l[k]) {
        i.first.l[k] = 0.;
        i.first.c[k]++;
      } else {
        i.second.l[k] = 1.;
        i.second.c[k]--;
      }
    } else if (i.second.c[k] < i.first.c[k]) {
      if (i.second.c[k] + 2 == i.first.c[k]) {
        i.second.l[k] = 0.;
        i.second.c[k]++;
        i.first.l[k] = 1.;
        i.first.c[k]--;
      } else if ((1. - i.second.l[k]) < i.first.l[k]) {
        i.second.l[k] = 0.;
        i.second.c[k]++;
      } else {
        i.first.l[k] = 1.;
        i.first.c[k]--;
      }
    }

    assert(i.first.c[k] == i.second.c[k]);
  }
}

template <int deg>
Scalar InterpolatedLevelSet<deg>::getTrilerp(const Coord &c) const {
  return m_phi.trilerp(c.c[0], c.c[1], c.c[2], c.l[0], c.l[1], c.l[2]);
}

template <int deg>
Scalar InterpolatedLevelSet<deg>::getInterpolated(const Coord &c) const {
  return getInterpolated(c.c[0], c.c[1], c.c[2], c.l[0], c.l[1], c.l[2]);
}

template <>
Scalar InterpolatedLevelSet<1>::getInterpolated(int i, int j, int k, Scalar fi,
                                                Scalar fj, Scalar fk) const {
  assert(i >= 0);
  assert(i < m_phi.ni);
  assert(j >= 0);
  assert(j < m_phi.nj);
  assert(k >= 0);
  assert(k < m_phi.nk);

  return m_phi.trilerp(i, j, k, fi, fj, fk);
}

template <>
Scalar InterpolatedLevelSet<3>::getInterpolated(int i, int j, int k, Scalar fi,
                                                Scalar fj, Scalar fk) const {
  assert(i >= 1);
  assert(i < m_phi.ni - 1);
  assert(j >= 1);
  assert(j < m_phi.nj - 1);
  assert(k >= 1);
  assert(k < m_phi.nk - 1);

  std::vector<Scalar> s(4), t(4), u(4);
  for (int di = -1; di <= 2; ++di) {
    for (int dj = -1; dj <= 2; ++dj) {
      for (int dk = -1; dk <= 2; ++dk) {
        s[dk + 1] = m_phi(i + di, j + dj, k + dk);
      }
      t[dj + 1] = bridson::cubic_interp(s[0], s[1], s[2], s[3], fk);
    }
    u[di + 1] = bridson::cubic_interp(t[0], t[1], t[2], t[3], fj);
  }
  return bridson::cubic_interp(u[0], u[1], u[2], u[3], fi);
}

template <>
Scalar InterpolatedLevelSet<1>::FractionInside(Scalar phi_left,
                                               Scalar phi_right) const {
  if (phi_left < 0 && phi_right < 0) return 1;
  if (phi_left < 0 && phi_right >= 0) return phi_left / (phi_left - phi_right);
  if (phi_left >= 0 && phi_right < 0)
    return phi_right / (phi_right - phi_left);
  else
    return 0;
}

template <>
Scalar InterpolatedLevelSet<3>::FractionInside(Scalar phi_left,
                                               Scalar phi_right) const {
  std::cerr << "Not implemented" << std::endl;
  return -1;
}

template <>
Scalar InterpolatedLevelSet<1>::FractionInside(Scalar phi_bl, Scalar phi_br,
                                               Scalar phi_tl,
                                               Scalar phi_tr) const {
  int inside_count = (phi_bl < 0 ? 1 : 0) + (phi_tl < 0 ? 1 : 0) +
                     (phi_br < 0 ? 1 : 0) + (phi_tr < 0 ? 1 : 0);
  Scalar list[] = {phi_bl, phi_br, phi_tr, phi_tl};

  if (inside_count == 4)
    return 1;
  else if (inside_count == 3) {
    // rotate until the positive value is in the first position
    while (list[0] < 0) {
      cycle_array(list, 4);
    }

    // Work out the area of the exterior triangle
    Scalar side0 = 1 - bridson::fraction_inside(list[0], list[3]);
    Scalar side1 = 1 - bridson::fraction_inside(list[0], list[1]);
    return 1 - 0.5 * side0 * side1;
  } else if (inside_count == 2) {
    // rotate until a negative value is in the first position, and the next
    // negative is in either slot 1 or 2.
    while (list[0] >= 0 || !(list[1] < 0 || list[2] < 0)) {
      cycle_array(list, 4);
    }

    if (list[1] < 0) {  // the matching signs are adjacent
      Scalar side_left = bridson::fraction_inside(list[0], list[3]);
      Scalar side_right = bridson::fraction_inside(list[1], list[2]);
      return 0.5 * (side_left + side_right);
    } else {  // matching signs are diagonally opposite
      // determine the centre point's sign to disambiguate this case
      Scalar middle_point = 0.25 * (list[0] + list[1] + list[2] + list[3]);
      if (middle_point < 0) {
        Scalar area = 0;

        // first triangle (top left)
        Scalar side1 = 1 - bridson::fraction_inside(list[0], list[3]);
        Scalar side3 = 1 - bridson::fraction_inside(list[2], list[3]);

        area += 0.5f * side1 * side3;

        // second triangle (top right)
        Scalar side2 = 1 - bridson::fraction_inside(list[2], list[1]);
        Scalar side0 = 1 - bridson::fraction_inside(list[0], list[1]);
        area += 0.5 * side0 * side2;

        return 1 - area;
      } else {
        Scalar area = 0;

        // first triangle (bottom left)
        Scalar side0 = bridson::fraction_inside(list[0], list[1]);
        Scalar side1 = bridson::fraction_inside(list[0], list[3]);
        area += 0.5 * side0 * side1;

        // second triangle (top right)
        Scalar side2 = bridson::fraction_inside(list[2], list[1]);
        Scalar side3 = bridson::fraction_inside(list[2], list[3]);
        area += 0.5 * side2 * side3;
        return area;
      }
    }
  } else if (inside_count == 1) {
    // rotate until the negative value is in the first position
    while (list[0] >= 0) {
      cycle_array(list, 4);
    }

    // Work out the area of the interior triangle, and subtract from 1.
    Scalar side0 = bridson::fraction_inside(list[0], list[3]);
    Scalar side1 = bridson::fraction_inside(list[0], list[1]);
    return 0.5 * side0 * side1;
  } else
    return 0;
}

template <>
Scalar InterpolatedLevelSet<3>::FractionInside(Scalar phi_bl, Scalar phi_br,
                                               Scalar phi_tl,
                                               Scalar phi_tr) const {
  std::cerr << "Not implemented" << std::endl;
  return -1;
}

template <int deg>
Scalar InterpolatedLevelSet<deg>::getFractionInsideU(const Vec3x &x) const {
  assert(m_initialized);

  // Compute the evaluation point in the level set's original coordinates
  const Vec4x &samplePoint =
      m_relativeTransformMatrix * Vec4x(x[0], x[1], x[2], 1.0f);

  // Compute the grid coordinates, integer and floating part
  int i, j, k;

  Scalar fi = (samplePoint[0] - m_origin[0]) / m_dx;
  Scalar fj = (samplePoint[1] - m_origin[1]) / m_dx;
  Scalar fk = (samplePoint[2] - m_origin[2]) / m_dx;

  bridson::get_barycentric(fi, i, fi, 1, m_phi.ni - 1);
  bridson::get_barycentric(fj, j, fj, 1, m_phi.nj - 1);
  bridson::get_barycentric(fk, k, fk, 1, m_phi.nk - 1);

  return getFractionInsideU(i, j, k, fi, fj, fk);
}

template <int deg>
Scalar InterpolatedLevelSet<deg>::getFractionInsideV(const Vec3x &x) const {
  assert(m_initialized);

  // Compute the evaluation point in the level set's original coordinates
  const Vec4x &samplePoint =
      m_relativeTransformMatrix * Vec4x(x[0], x[1], x[2], 1.0f);

  // Compute the grid coordinates, integer and floating part
  int i, j, k;

  Scalar fi = (samplePoint[0] - m_origin[0]) / m_dx;
  Scalar fj = (samplePoint[1] - m_origin[1]) / m_dx;
  Scalar fk = (samplePoint[2] - m_origin[2]) / m_dx;

  bridson::get_barycentric(fi, i, fi, 1, m_phi.ni - 1);
  bridson::get_barycentric(fj, j, fj, 1, m_phi.nj - 1);
  bridson::get_barycentric(fk, k, fk, 1, m_phi.nk - 1);

  return getFractionInsideV(i, j, k, fi, fj, fk);
}

template <int deg>
Scalar InterpolatedLevelSet<deg>::getFractionInsideW(const Vec3x &x) const {
  assert(m_initialized);

  // Compute the evaluation point in the level set's original coordinates
  const Vec4x &samplePoint =
      m_relativeTransformMatrix * Vec4x(x[0], x[1], x[2], 1.0f);

  // Compute the grid coordinates, integer and floating part
  int i, j, k;

  Scalar fi = (samplePoint[0] - m_origin[0]) / m_dx;
  Scalar fj = (samplePoint[1] - m_origin[1]) / m_dx;
  Scalar fk = (samplePoint[2] - m_origin[2]) / m_dx;

  bridson::get_barycentric(fi, i, fi, 1, m_phi.ni - 1);
  bridson::get_barycentric(fj, j, fj, 1, m_phi.nj - 1);
  bridson::get_barycentric(fk, k, fk, 1, m_phi.nk - 1);

  return getFractionInsideW(i, j, k, fi, fj, fk);
}

template <>
Scalar InterpolatedLevelSet<1>::getFractionInsideU(int i, int j, int k,
                                                   Scalar fi, Scalar fj,
                                                   Scalar fk) const {
  Scalar phi00 = getInterpolated(i, j, k, fi, 0.0f, 0.0f);
  Scalar phi10 = getInterpolated(i, j, k, fi, 1.0f, 0.0f);
  Scalar phi01 = getInterpolated(i, j, k, fi, 0.0f, 1.0f);
  Scalar phi11 = getInterpolated(i, j, k, fi, 1.0f, 1.0f);

  return FractionInside(phi00, phi10, phi01, phi11);
}

template <>
Scalar InterpolatedLevelSet<1>::getFractionInsideV(int i, int j, int k,
                                                   Scalar fi, Scalar fj,
                                                   Scalar fk) const {
  Scalar phi00 = getInterpolated(i, j, k, 0.0f, fj, 0.0f);
  Scalar phi10 = getInterpolated(i, j, k, 1.0f, fj, 0.0f);
  Scalar phi01 = getInterpolated(i, j, k, 0.0f, fj, 1.0f);
  Scalar phi11 = getInterpolated(i, j, k, 1.0f, fj, 1.0f);

  return FractionInside(phi00, phi10, phi01, phi11);
}

template <>
Scalar InterpolatedLevelSet<1>::getFractionInsideW(int i, int j, int k,
                                                   Scalar fi, Scalar fj,
                                                   Scalar fk) const {
  Scalar phi00 = getInterpolated(i, j, k, 0.0f, 0.0f, fk);
  Scalar phi10 = getInterpolated(i, j, k, 1.0f, 0.0f, fk);
  Scalar phi01 = getInterpolated(i, j, k, 0.0f, 1.0f, fk);
  Scalar phi11 = getInterpolated(i, j, k, 1.0f, 1.0f, fk);

  return FractionInside(phi00, phi10, phi01, phi11);
}

template <>
Scalar InterpolatedLevelSet<3>::getFractionInsideU(int i, int j, int k,
                                                   Scalar fi, Scalar fj,
                                                   Scalar fk) const {
  std::cerr << "Not implemented" << std::endl;
  return -1;
}

template <>
Scalar InterpolatedLevelSet<3>::getFractionInsideV(int i, int j, int k,
                                                   Scalar fi, Scalar fj,
                                                   Scalar fk) const {
  std::cerr << "Not implemented" << std::endl;
  return -1;
}

template <>
Scalar InterpolatedLevelSet<3>::getFractionInsideW(int i, int j, int k,
                                                   Scalar fi, Scalar fj,
                                                   Scalar fk) const {
  std::cerr << "Not implemented" << std::endl;
  return -1;
}

template <>
Vec3x InterpolatedLevelSet<1>::getInterpolatedGradient(int i, int j, int k,
                                                       Scalar fi, Scalar fj,
                                                       Scalar fk) const {
  Vec3x phi0, phi1;

  phi0[0] = getInterpolated(i, j, k, 0.0f, fj, fk);
  phi1[0] = getInterpolated(i, j, k, 1.0f, fj, fk);
  phi0[1] = getInterpolated(i, j, k, fi, 0.0f, fk);
  phi1[1] = getInterpolated(i, j, k, fi, 1.0f, fk);
  phi0[2] = getInterpolated(i, j, k, fi, fj, 0.0f);
  phi1[2] = getInterpolated(i, j, k, fi, fj, 1.0f);

  return (phi1 - phi0) / m_dx;
}

template <>
Vec3x InterpolatedLevelSet<3>::getInterpolatedGradient(int i, int j, int k,
                                                       Scalar fi, Scalar fj,
                                                       Scalar fk) const {
  assert(i >= 1);
  assert(i < m_phi.ni - 1);
  assert(j >= 1);
  assert(j < m_phi.nj - 1);
  assert(k >= 1);
  assert(k < m_phi.nk - 1);

  Vec3x gradf;
  std::vector<Scalar> s(4), t(4), u(4);

  for (int di = -1; di <= 2; ++di) {
    for (int dj = -1; dj <= 2; ++dj) {
      for (int dk = -1; dk <= 2; ++dk) {
        s[dk + 1] = m_phi(i + di, j + dj, k + dk);
      }
      t[dj + 1] = bridson::cubic_interp_diff(s[0], s[1], s[2], s[3], fk);
    }
    u[di + 1] = bridson::cubic_interp(t[0], t[1], t[2], t[3], fj);
  }
  gradf[2] = bridson::cubic_interp(u[0], u[1], u[2], u[3], fi);

  for (int di = -1; di <= 2; ++di) {
    for (int dj = -1; dj <= 2; ++dj) {
      for (int dk = -1; dk <= 2; ++dk) {
        s[dk + 1] = m_phi(i + di, j + dj, k + dk);
      }
      t[dj + 1] = bridson::cubic_interp(s[0], s[1], s[2], s[3], fk);
    }
    u[di + 1] = bridson::cubic_interp_diff(t[0], t[1], t[2], t[3], fj);
  }
  gradf[1] = bridson::cubic_interp(u[0], u[1], u[2], u[3], fi);

  for (int di = -1; di <= 2; ++di) {
    for (int dj = -1; dj <= 2; ++dj) {
      for (int dk = -1; dk <= 2; ++dk) {
        s[dk + 1] = m_phi(i + di, j + dj, k + dk);
      }
      t[dj + 1] = bridson::cubic_interp(s[0], s[1], s[2], s[3], fk);
    }
    u[di + 1] = bridson::cubic_interp(t[0], t[1], t[2], t[3], fj);
  }
  gradf[0] = bridson::cubic_interp_diff(u[0], u[1], u[2], u[3], fi);

  return gradf / m_dx;
}

template <>
Mat3x InterpolatedLevelSet<1>::getInterpolatedHessian(int i, int j, int k,
                                                      Scalar fi, Scalar fj,
                                                      Scalar fk) const {
  assert(i >= 1);
  assert(i < m_phi.ni - 1);
  assert(j >= 1);
  assert(j < m_phi.nj - 1);
  assert(k >= 1);
  assert(k < m_phi.nk - 1);

  Mat3x hessf;
  hessf.setZero();

  hessf(1, 0) = hessf(0, 1) =
      getInterpolated(i, j, k, 1, 1, fk) - getInterpolated(i, j, k, 1, 0, fk) -
      getInterpolated(i, j, k, 0, 1, fk) + getInterpolated(i, j, k, 0, 0, fk);
  hessf(2, 0) = hessf(0, 2) =
      getInterpolated(i, j, k, 1, fj, 1) - getInterpolated(i, j, k, 1, fj, 0) -
      getInterpolated(i, j, k, 0, fj, 1) + getInterpolated(i, j, k, 0, fj, 0);
  hessf(2, 1) = hessf(1, 2) =
      getInterpolated(i, j, k, fi, 1, 1) - getInterpolated(i, j, k, fi, 1, 0) -
      getInterpolated(i, j, k, fi, 0, 1) + getInterpolated(i, j, k, fi, 0, 0);

  return hessf / (m_dx * m_dx);
}

template <>
Mat3x InterpolatedLevelSet<3>::getInterpolatedHessian(int i, int j, int k,
                                                      Scalar fi, Scalar fj,
                                                      Scalar fk) const {
  assert(i >= 1);
  assert(i < m_phi.ni - 1);
  assert(j >= 1);
  assert(j < m_phi.nj - 1);
  assert(k >= 1);
  assert(k < m_phi.nk - 1);

  Mat3x hessf;
  std::vector<Scalar> s(4), t(4), u(4);

  for (int di = -1; di <= 2; ++di) {
    for (int dj = -1; dj <= 2; ++dj) {
      for (int dk = -1; dk <= 2; ++dk) {
        s[dk + 1] = m_phi(i + di, j + dj, k + dk);
      }
      t[dj + 1] = bridson::cubic_interp_diff2(s[0], s[1], s[2], s[3], fk);
    }
    u[di + 1] = bridson::cubic_interp(t[0], t[1], t[2], t[3], fj);
  }
  hessf(2, 2) = bridson::cubic_interp(u[0], u[1], u[2], u[3], fi);

  for (int di = -1; di <= 2; ++di) {
    for (int dj = -1; dj <= 2; ++dj) {
      for (int dk = -1; dk <= 2; ++dk) {
        s[dk + 1] = m_phi(i + di, j + dj, k + dk);
      }
      t[dj + 1] = bridson::cubic_interp(s[0], s[1], s[2], s[3], fk);
    }
    u[di + 1] = bridson::cubic_interp_diff2(t[0], t[1], t[2], t[3], fj);
  }
  hessf(1, 1) = bridson::cubic_interp(u[0], u[1], u[2], u[3], fi);

  for (int di = -1; di <= 2; ++di) {
    for (int dj = -1; dj <= 2; ++dj) {
      for (int dk = -1; dk <= 2; ++dk) {
        s[dk + 1] = m_phi(i + di, j + dj, k + dk);
      }
      t[dj + 1] = bridson::cubic_interp(s[0], s[1], s[2], s[3], fk);
    }
    u[di + 1] = bridson::cubic_interp(t[0], t[1], t[2], t[3], fj);
  }
  hessf(0, 0) = bridson::cubic_interp_diff2(u[0], u[1], u[2], u[3], fi);

  for (int di = -1; di <= 2; ++di) {
    for (int dj = -1; dj <= 2; ++dj) {
      for (int dk = -1; dk <= 2; ++dk) {
        s[dk + 1] = m_phi(i + di, j + dj, k + dk);
      }
      t[dj + 1] = bridson::cubic_interp_diff(s[0], s[1], s[2], s[3], fk);
    }
    u[di + 1] = bridson::cubic_interp(t[0], t[1], t[2], t[3], fj);
  }
  hessf(0, 2) = hessf(2, 0) =
      bridson::cubic_interp_diff(u[0], u[1], u[2], u[3], fi);

  for (int di = -1; di <= 2; ++di) {
    for (int dj = -1; dj <= 2; ++dj) {
      for (int dk = -1; dk <= 2; ++dk) {
        s[dk + 1] = m_phi(i + di, j + dj, k + dk);
      }
      t[dj + 1] = bridson::cubic_interp_diff(s[0], s[1], s[2], s[3], fk);
    }
    u[di + 1] = bridson::cubic_interp_diff(t[0], t[1], t[2], t[3], fj);
  }
  hessf(1, 2) = hessf(2, 1) = bridson::cubic_interp(u[0], u[1], u[2], u[3], fi);

  for (int di = -1; di <= 2; ++di) {
    for (int dj = -1; dj <= 2; ++dj) {
      for (int dk = -1; dk <= 2; ++dk) {
        s[dk + 1] = m_phi(i + di, j + dj, k + dk);
      }
      t[dj + 1] = bridson::cubic_interp(s[0], s[1], s[2], s[3], fk);
    }
    u[di + 1] = bridson::cubic_interp_diff(t[0], t[1], t[2], t[3], fj);
  }
  hessf(0, 1) = hessf(1, 0) =
      bridson::cubic_interp_diff(u[0], u[1], u[2], u[3], fi);

  return hessf / (m_dx * m_dx);
}

void computeLevelSet(const FaceArray &triangles,
                     const std::vector<unsigned> &triIndices,
                     const Vec3xArray &x, const Vec3xArray &v, Vec3x &origin,
                     const Vec3x &center, Scalar &dx, bridson::Array3x &phi,
                     bridson::Array3<Vec3x, bridson::Array1<Vec3x> > &phiVel,
                     bridson::Array3i &closest_tri, bool adaptativeDx) {
  const int exact_band = 1;

  const Vec3x exactOrigin = origin;  // origin will be padded later
  const Vec3x length = 2 * (center - exactOrigin);
  const Scalar minDxCube =
      length[0] * length[1] * length[2] / MAX_LEVEL_SET_CELLS;

  bool acceptDx = true;
  int ni, nj, nk;

  bridson::Array3i i_in_intersection_count(
      1, 1, 1, 0);  // intersection_count(i,j,k) is # of tri intersections in
                    // (i-1,i]x{j}x{k}
  bridson::Array3i i_out_intersection_count(
      1, 1, 1, 0);  // intersection_count(i,j,k) is # of tri intersections in
                    // (i-1,i]x{j}x{k}
  bridson::Array3i j_in_intersection_count(
      1, 1, 1, 0);  // intersection_count(i,j,k) is # of tri intersections in
                    // (i-1,i]x{j}x{k}
  bridson::Array3i j_out_intersection_count(
      1, 1, 1, 0);  // intersection_count(i,j,k) is # of tri intersections in
                    // (i-1,i]x{j}x{k}
  bridson::Array3i k_in_intersection_count(
      1, 1, 1, 0);  // intersection_count(i,j,k) is # of tri intersections in
                    // (i-1,i]x{j}x{k}
  bridson::Array3i k_out_intersection_count(
      1, 1, 1, 0);  // intersection_count(i,j,k) is # of tri intersections in
                    // (i-1,i]x{j}x{k}

  do {
    // Padding (to get nice zeros on border)
    origin = exactOrigin - 2.5 * dx * Vec3x::Ones();

    ni = std::ceil(length[0] / dx) + 6;
    nj = std::ceil(length[1] / dx) + 6;
    nk = std::ceil(length[2] / dx) + 6;

    phi.resize(ni, nj, nk);
    phiVel.resize(ni, nj, nk);
    phi.assign((ni + nj + nk) * dx);  // upper bound on distance

    i_in_intersection_count.resize(ni, nj, nk);
    i_out_intersection_count.resize(ni, nj, nk);
    j_in_intersection_count.resize(ni, nj, nk);
    j_out_intersection_count.resize(ni, nj, nk);
    k_in_intersection_count.resize(ni, nj, nk);
    k_out_intersection_count.resize(ni, nj, nk);
    i_in_intersection_count.assign(0);
    i_out_intersection_count.assign(0);
    j_in_intersection_count.assign(0);
    j_out_intersection_count.assign(0);
    k_in_intersection_count.assign(0);
    k_out_intersection_count.assign(0);

    closest_tri.resize(ni, nj, nk);
    closest_tri.assign(-1);

    // we begin by initializing distances near the mesh, and figuring out
    // intersection counts
    //
    Vec3x ijkmin, ijkmax;

    //    for(unsigned int t=0; t<triangles.size(); ++t){
    for (std::vector<unsigned>::const_iterator tItr = triIndices.begin();
         tItr != triIndices.end(); ++tItr) {
      unsigned t = *tItr;
      //        unsigned int p, q, r; assign(triangles[t], p, q, r);
      const unsigned int p = triangles[t](0), q = triangles[t](1),
                         r = triangles[t](2);

      // coordinates in grid to high precision
      const Scalar fip = ((Scalar)x[p][0] - origin[0]) / dx,
                   fjp = ((Scalar)x[p][1] - origin[1]) / dx,
                   fkp = ((Scalar)x[p][2] - origin[2]) / dx;
      const Scalar fiq = ((Scalar)x[q][0] - origin[0]) / dx,
                   fjq = ((Scalar)x[q][1] - origin[1]) / dx,
                   fkq = ((Scalar)x[q][2] - origin[2]) / dx;
      const Scalar fir = ((Scalar)x[r][0] - origin[0]) / dx,
                   fjr = ((Scalar)x[r][1] - origin[1]) / dx,
                   fkr = ((Scalar)x[r][2] - origin[2]) / dx;

      //        // do distances nearby
      int i0 = clamp(int(bridson::min3(fip, fiq, fir)) - exact_band, 0, ni - 1),
          i1 = clamp(int(bridson::max3(fip, fiq, fir)) + exact_band + 1, 0,
                     ni - 1);
      int j0 = clamp(int(bridson::min3(fjp, fjq, fjr)) - exact_band, 0, nj - 1),
          j1 = clamp(int(bridson::max3(fjp, fjq, fjr)) + exact_band + 1, 0,
                     nj - 1);
      int k0 = clamp(int(bridson::min3(fkp, fkq, fkr)) - exact_band, 0, nk - 1),
          k1 = clamp(int(bridson::max3(fkp, fkq, fkr)) + exact_band + 1, 0,
                     nk - 1);

      for (int k = k0; k <= k1; ++k)
        for (int j = j0; j <= j1; ++j)
          for (int i = i0; i <= i1; ++i) {
            Vec3x gx(i * dx + origin[0], j * dx + origin[1],
                     k * dx + origin[2]);
            Scalar t1, t2, t3;
            Scalar d =
                point_triangle_distance(gx, x[p], x[q], x[r], t1, t2, t3);
            if (d < phi(i, j, k)) {
              phi(i, j, k) = d;
              phiVel(i, j, k) = (v[p] * t1 + v[q] * t2 + v[r] * t3);
              //                        convert3( phiVel( i, j, k ), ( v[p] * t1
              //                        + v[q] * t2 + v[r] * t3 ) );
              closest_tri(i, j, k) = t;
            }
          }

      // and do intersection counts
      i0 = clamp((int)std::ceil(bridson::min3(fip, fiq, fir)), 0, ni - 1);
      i1 = clamp((int)std::floor(bridson::max3(fip, fiq, fir)), 0, ni - 1);
      j0 = clamp((int)std::ceil(bridson::min3(fjp, fjq, fjr)), 0, nj - 1);
      j1 = clamp((int)std::floor(bridson::max3(fjp, fjq, fjr)), 0, nj - 1);
      k0 = clamp((int)std::ceil(bridson::min3(fkp, fkq, fkr)), 0, nk - 1);
      k1 = clamp((int)std::floor(bridson::max3(fkp, fkq, fkr)), 0, nk - 1);

      // Ok, this is really repetitive and could be factored but, but not easily
      // because of this whole
      // i/j/k naming convention

      const Scalar fin =
          ((fjq - fjp) * (fkr - fkp)) - ((fjr - fjp) * (fkq - fkp));
      {
        bridson::Array3i &intersection_count =
            fin > 0. ? i_out_intersection_count : i_in_intersection_count;

        for (int k = k0; k <= k1; ++k)
          for (int j = j0; j <= j1; ++j) {
            Scalar a, b, c;
            if (point_in_triangle_2d(j, k, fjp, fkp, fjq, fkq, fjr, fkr, a, b,
                                     c)) {
              Scalar fi =
                  a * fip + b * fiq + c * fir;  // intersection i coordinate
              int i_interval = int(std::ceil(
                  fi));  // intersection is in (i_interval-1,i_interval]
              if (i_interval < 0)
                ++intersection_count(
                    0, j, k);  // we enlarge the first interval to include
                               // everything to the -x direction

              else if (i_interval < ni)
                ++intersection_count(i_interval, j, k);
              // we ignore intersections that are beyond the +x side of the grid
            }
          }
      }

      const Scalar fjn =
          ((fkq - fkp) * (fir - fip)) - ((fkr - fkp) * (fiq - fip));
      {
        bridson::Array3i &intersection_count =
            fjn > 0. ? j_out_intersection_count : j_in_intersection_count;
        for (int i = i0; i <= i1; ++i)
          for (int k = k0; k <= k1; ++k) {
            Scalar a, b, c;
            if (point_in_triangle_2d(k, i, fkp, fip, fkq, fiq, fkr, fir, a, b,
                                     c)) {
              Scalar fj =
                  a * fjp + b * fjq + c * fjr;  // intersection i coordinate
              int j_interval = int(std::ceil(
                  fj));  // intersection is in (i_interval-1,i_interval]
              if (j_interval < 0)
                ++intersection_count(
                    i, 0, k);  // we enlarge the first interval to include
                               // everything to the -x direction

              else if (j_interval < nj)
                ++intersection_count(i, j_interval, k);
              // we ignore intersections that are beyond the +x side of the grid
            }
          }
      }

      const Scalar fkn =
          -((fjq - fjp) * (fir - fip)) + ((fjr - fjp) * (fiq - fip));
      {
        bridson::Array3i &intersection_count =
            fkn > 0. ? k_out_intersection_count : k_in_intersection_count;
        for (int i = i0; i <= i1; ++i)
          for (int j = j0; j <= j1; ++j) {
            Scalar a, b, c;
            if (point_in_triangle_2d(j, i, fjp, fip, fjq, fiq, fjr, fir, a, b,
                                     c)) {
              Scalar fk =
                  a * fkp + b * fkq + c * fkr;  // intersection i coordinate
              int k_interval = int(std::ceil(
                  fk));  // intersection is in (i_interval-1,i_interval]
              if (k_interval < 0)
                ++intersection_count(
                    i, j, 0);  // we enlarge the first interval to include
                               // everything to the -x direction

              else if (k_interval < nk)
                ++intersection_count(i, j, k_interval);
              // we ignore intersections that are beyond the +x side of the grid
            }
          }
      }
    }

    if (adaptativeDx) {
      // Interaction count stats
      std::vector<unsigned> histogram;

      unsigned totInts = 0;
      unsigned intersectingCells = 0;
      for (int i = 0; i < ni; ++i) {
        for (int j = 0; j < nj; ++j) {
          for (int k = 0; k < nk; ++k) {
            unsigned totCellInt =
                std::max(i_in_intersection_count(i, j, k) +
                             i_out_intersection_count(i, j, k),
                         std::max(j_in_intersection_count(i, j, k) +
                                      j_out_intersection_count(i, j, k),
                                  k_in_intersection_count(i, j, k) +
                                      k_out_intersection_count(i, j, k)));

            totInts += totCellInt;

            if (totCellInt) {
              if (histogram.size() < totCellInt) {
                histogram.resize(totCellInt, 0);
              }
              ++histogram[totCellInt - 1];
              ++intersectingCells;
            }
          }
        }
      }

      if (totInts) {
        const Scalar propOfUnique =
            (((Scalar)histogram[0]) / intersectingCells);
        //                std::cout << dx << " -> Proportion of unique
        //                intersects: " <<  propOfUnique << std::endl ;

        acceptDx = (propOfUnique > MIN_UNIQUE_INTERSECTIONS_RATIO);

        if (!acceptDx) {
          const Scalar newDx = dx * .5f;
          acceptDx = (newDx * newDx * newDx < minDxCube);

          if (!acceptDx) dx = newDx;
        }
      }
    }

  } while (!acceptDx);

  // FIXME find deterministic parallelization scheme

  // and now we fill in the rest of the distances with fast sweeping
  //    omp_set_num_threads( 4 );
  for (unsigned int pass = 0; pass < 2; ++pass) {
#pragma omp parallel sections
    {
#pragma omp section
      {
        sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, +1, +1,
              +1);
      }
#pragma omp section
      {
        sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, -1, -1,
              -1);
      }
#pragma omp section
      {
        sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, +1, +1,
              -1);
      }
#pragma omp section
      {
        sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, -1, -1,
              +1);
      }
#pragma omp section
      {
        sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, +1, -1,
              +1);
      }
#pragma omp section
      {
        sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, -1, +1,
              -1);
      }
#pragma omp section
      {
        sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, +1, -1,
              -1);
      }
#pragma omp section
      {
        sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, -1, +1,
              +1);
      }
    }
  }

  bridson::Array3i interior_votes(ni, nj, nk, 0);

  // then figure out signs (inside/outside) from intersection counts
  // (see above)
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < nj; ++j) {
      std::vector<int> total_count(ni, 0);

      int cur_count = 0;
      int min_count = 0;

      for (int i = 0; i < ni; ++i) {
        cur_count += i_in_intersection_count(i, j, k) -
                     i_out_intersection_count(i, j, k);

        if (cur_count < min_count) {
          min_count = cur_count;
        }

        total_count[i] = cur_count;
      }
      for (int i = 0; i < ni; ++i) {
        if (total_count[i] > min_count) {
          ++interior_votes(i, j, k);
        }
      }
    }

    for (int i = 0; i < ni; ++i) {
      std::vector<int> total_count(nj, 0);

      int cur_count = 0;
      int min_count = 0;

      for (int j = 0; j < nj; ++j) {
        cur_count += j_in_intersection_count(i, j, k) -
                     j_out_intersection_count(i, j, k);

        if (cur_count < min_count) {
          min_count = cur_count;
        }

        total_count[j] = cur_count;
      }
      for (int j = 0; j < nj; ++j) {
        if (total_count[j] > min_count) {
          ++interior_votes(i, j, k);
        }
      }
    }
  }

  // then figure out signs (inside/outside) from intersection counts
  for (int i = 0; i < ni; ++i) {
    for (int j = 0; j < nj; ++j) {
      std::vector<int> total_count(nk, 0);

      int cur_count = 0;
      int min_count = 0;

      for (int k = 0; k < nk; ++k) {
        cur_count += k_in_intersection_count(i, j, k) -
                     k_out_intersection_count(i, j, k);

        if (cur_count < min_count) {
          min_count = cur_count;
        }

        total_count[k] = cur_count;
      }
      for (int k = 0; k < nk; ++k) {
        if (total_count[k] > min_count) {
          ++interior_votes(i, j, k);
        }
      }
    }
  }

  // If at least two directions agree on an interior point, change level set
  // sign
  for (int i = 0; i < ni; ++i)
    for (int j = 0; j < nj; ++j) {
      for (int k = 0; k < nk; ++k) {
        if (interior_votes(i, j, k) > 1) {
          phi(i, j, k) *= -1;
        }
      }
    }
}

void check_neighbour(const FaceArray &tri, const Vec3xArray &x,
                     const Vec3xArray &v, bridson::Array3x &phi,
                     bridson::Array3<Vec3x, bridson::Array1<Vec3x> > &phi_vel,
                     bridson::Array3i &closest_tri, const Vec3x &gx,
                     const int i0, const int j0, const int k0, const int i1,
                     const int j1, const int k1) {
  if (closest_tri(i1, j1, k1) >= 0) {
    //        unsigned int p, q, r; assign(tri[closest_tri(i1,j1,k1)], p, q, r);
    unsigned int p = tri[closest_tri(i1, j1, k1)](0);
    unsigned int q = tri[closest_tri(i1, j1, k1)](1);
    unsigned int r = tri[closest_tri(i1, j1, k1)](2);

    Scalar t1, t2, t3;
    Scalar d = point_triangle_distance(gx, x[p], x[q], x[r], t1, t2, t3);
    if (d < phi(i0, j0, k0)) {
      phi(i0, j0, k0) = d;
      phi_vel(i0, j0, k0) = (v[p] * t1 + v[q] * t2 + v[r] * t3);
      //            convert3( phi_vel( i0, j0, k0 ), ( v[p] * t1 + v[q] * t2 +
      //            v[r] * t3 ) );
      closest_tri(i0, j0, k0) = closest_tri(i1, j1, k1);
    }
  }
}

void sweep(const FaceArray &tri, const Vec3xArray &x, const Vec3xArray &v,
           bridson::Array3x &phi,
           bridson::Array3<Vec3x, bridson::Array1<Vec3x> > &phi_vel,
           bridson::Array3i &closest_tri, const Vec3x &origin, const Scalar dx,
           const int di, const int dj, const int dk) {
  int i0, i1;
  if (di > 0) {
    i0 = 1;
    i1 = phi.ni;
  } else {
    i0 = phi.ni - 2;
    i1 = -1;
  }
  int j0, j1;
  if (dj > 0) {
    j0 = 1;
    j1 = phi.nj;
  } else {
    j0 = phi.nj - 2;
    j1 = -1;
  }
  int k0, k1;
  if (dk > 0) {
    k0 = 1;
    k1 = phi.nk;
  } else {
    k0 = phi.nk - 2;
    k1 = -1;
  }
  for (int k = k0; k != k1; k += dk)
    for (int j = j0; j != j1; j += dj)
      for (int i = i0; i != i1; i += di) {
        Vec3x gx(i * dx + origin[0], j * dx + origin[1], k * dx + origin[2]);
        check_neighbour(tri, x, v, phi, phi_vel, closest_tri, gx, i, j, k,
                        i - di, j, k);
        check_neighbour(tri, x, v, phi, phi_vel, closest_tri, gx, i, j, k,
                        i - di, j - dj, k);
        check_neighbour(tri, x, v, phi, phi_vel, closest_tri, gx, i, j, k, i, j,
                        k - dk);
        check_neighbour(tri, x, v, phi, phi_vel, closest_tri, gx, i, j, k,
                        i - di, j, k - dk);
        check_neighbour(tri, x, v, phi, phi_vel, closest_tri, gx, i, j, k, i,
                        j - dj, k - dk);
        check_neighbour(tri, x, v, phi, phi_vel, closest_tri, gx, i, j, k,
                        i - di, j - dj, k - dk);
      }
}

int orientation(const Scalar x1, const Scalar y1, const Scalar x2,
                const Scalar y2, Scalar &twice_signed_area) {
  twice_signed_area = y1 * x2 - x1 * y2;
  if (twice_signed_area > 0)
    return 1;
  else if (twice_signed_area < 0)
    return -1;
  else if (y2 > y1)
    return 1;
  else if (y2 < y1)
    return -1;
  else if (x1 > x2)
    return 1;
  else if (x1 < x2)
    return -1;
  else
    return 0;  // only true when x1==x2 and y1==y2
}

// Checks whether the point (x0, y0) is in the triangle; if so sets a, b, c to
// barycentric coordinates
bool point_in_triangle_2d(const Scalar x0, const Scalar y0, Scalar x1,
                          Scalar y1, Scalar x2, Scalar y2, Scalar x3, Scalar y3,
                          Scalar &a, Scalar &b, Scalar &c) {
  x1 -= x0;
  x2 -= x0;
  x3 -= x0;
  y1 -= y0;
  y2 -= y0;
  y3 -= y0;

  int signa = orientation(x2, y2, x3, y3, a);
  if (signa == 0) return false;

  int signb = orientation(x3, y3, x1, y1, b);
  if (signb != signa) return false;

  int signc = orientation(x1, y1, x2, y2, c);
  if (signc != signa) return false;

  Scalar sum = a + b + c;
  assert(sum != 0);
  // if the SOS signs match and are nonzero, there's no way all of a, b, and c
  // are zero.
  a /= sum;
  b /= sum;
  c /= sum;
  return true;
}

void computeTrilinearValueCoeffs(const bridson::Array3x &phi,
                                 const InterpolatedLevelSet<1>::Interval &iv,
                                 Vec4x &coeffs) {
  coeffs.setZero();

  Vec3x lambdaA[2];
  Vec3x lambdaB[2];

  // lambda(t) = lambdaA[1] t + lambdaB[1]
  lambdaA[1] = (iv.second.l - iv.first.l) / (iv.second.t - iv.first.t);
  lambdaB[1] = iv.first.l - iv.first.t * lambdaA[1];
  // (1 - lambda(t)) = lambdaA[0] t + lambdaB[0]
  lambdaA[0] = -lambdaA[1];
  lambdaB[0] = Vec3x::Ones() - lambdaB[1];

  for (unsigned i = 0; i < 2; ++i)
    for (unsigned j = 0; j < 2; ++j)
      for (unsigned k = 0; k < 2; ++k) {
        const Scalar val =
            phi(iv.first.c[0] + i, iv.first.c[1] + j, iv.first.c[2] + k);
        coeffs[3] += val * (lambdaA[i][0] * lambdaA[j][1] * lambdaA[k][2]);
        coeffs[2] += val * (lambdaA[i][0] * lambdaA[j][1] * lambdaB[k][2] +
                            lambdaB[i][0] * lambdaA[j][1] * lambdaA[k][2] +
                            lambdaA[i][0] * lambdaB[j][1] * lambdaA[k][2]);
        coeffs[1] += val * (lambdaA[i][0] * lambdaB[j][1] * lambdaB[k][2] +
                            lambdaB[i][0] * lambdaB[j][1] * lambdaA[k][2] +
                            lambdaB[i][0] * lambdaA[j][1] * lambdaB[k][2]);
        coeffs[0] += val * (lambdaB[i][0] * lambdaB[j][1] * lambdaB[k][2]);
      }
}

void computeBilinearGradientCoeffs(const bridson::Array3x &phi,
                                   const InterpolatedLevelSet<1>::Interval &iv,
                                   Mat3x &coeffs) {
  coeffs.setZero();

  Vec3x lambdaA[2];
  Vec3x lambdaB[2];

  // lambda(t) = lambdaA[1] t + lambdaB[1]
  lambdaA[1] = (iv.second.l - iv.first.l) / (iv.second.t - iv.first.t);
  lambdaB[1] = iv.first.l - iv.first.t * lambdaA[1];
  // (1 - lambda(t)) = lambdaA[0] t + lambdaB[0]
  lambdaA[0] = -lambdaA[1];
  lambdaB[0] = Vec3x::Ones() - lambdaB[1];

  for (unsigned j = 0; j < 2; ++j)
    for (unsigned k = 0; k < 2; ++k) {
      const Scalar grad =
          phi(iv.first.c[0] + 1, iv.first.c[1] + j, iv.first.c[2] + k) -
          phi(iv.first.c[0], iv.first.c[1] + j, iv.first.c[2] + k);

      coeffs(0, 2) += grad * (lambdaA[j][1] * lambdaA[k][2]);
      coeffs(0, 1) += grad * (lambdaA[j][1] * lambdaB[k][2] +
                              lambdaB[j][1] * lambdaA[k][2]);
      coeffs(0, 0) += grad * (lambdaB[j][1] * lambdaB[k][2]);
    }
  for (unsigned i = 0; i < 2; ++i)
    for (unsigned k = 0; k < 2; ++k) {
      const Scalar grad =
          phi(iv.first.c[0] + i, iv.first.c[1] + 1, iv.first.c[2] + k) -
          phi(iv.first.c[0] + i, iv.first.c[1], iv.first.c[2] + k);

      coeffs(1, 2) += grad * (lambdaA[i][0] * lambdaA[k][2]);
      coeffs(1, 1) += grad * (lambdaA[i][0] * lambdaB[k][2] +
                              lambdaB[i][0] * lambdaA[k][2]);
      coeffs(1, 0) += grad * (lambdaB[i][0] * lambdaB[k][2]);
    }
  for (unsigned j = 0; j < 2; ++j)
    for (unsigned i = 0; i < 2; ++i) {
      const Scalar grad =
          phi(iv.first.c[0] + i, iv.first.c[1] + j, iv.first.c[2] + 1) -
          phi(iv.first.c[0] + i, iv.first.c[1] + j, iv.first.c[2]);

      coeffs(2, 2) += grad * (lambdaA[j][1] * lambdaA[i][0]);
      coeffs(2, 1) += grad * (lambdaA[j][1] * lambdaB[i][0] +
                              lambdaB[j][1] * lambdaA[i][0]);
      coeffs(2, 0) += grad * (lambdaB[j][1] * lambdaB[i][0]);
    }
}

void computeLinearHessianCoeffs(const bridson::Array3x &phi,
                                const InterpolatedLevelSet<1>::Interval &iv,
                                Mat3x &HA, Mat3x &HB) {
  HA.setZero();
  HB.setZero();

  Vec3x lambdaA[2];
  Vec3x lambdaB[2];

  // lambda(t) = lambdaA[1] t + lambdaB[1]
  lambdaA[1] = (iv.second.l - iv.first.l) / (iv.second.t - iv.first.t);
  lambdaB[1] = iv.first.l - iv.first.t * lambdaA[1];
  // (1 - lambda(t)) = lambdaA[0] t + lambdaB[0]
  lambdaA[0] = -lambdaA[1];
  lambdaB[0] = Vec3x::Ones() - lambdaB[1];

  for (unsigned i = 0; i < 2; ++i) {
    const Scalar h =
        phi(iv.first.c[0] + i, iv.first.c[1] + 1, iv.first.c[2] + 1) -
        phi(iv.first.c[0] + i, iv.first.c[1] + 1, iv.first.c[2] + 0) -
        phi(iv.first.c[0] + i, iv.first.c[1] + 0, iv.first.c[2] + 1) +
        phi(iv.first.c[0] + i, iv.first.c[1] + 0, iv.first.c[2] + 0);
    HA(2, 1) += lambdaA[i][0] * h;
    HA(1, 2) += lambdaA[i][0] * h;
    HB(2, 1) += lambdaB[i][0] * h;
    HB(1, 2) += lambdaB[i][0] * h;
  }
  for (unsigned j = 0; j < 2; ++j) {
    const Scalar h =
        phi(iv.first.c[0] + 1, iv.first.c[1] + j, iv.first.c[2] + 1) -
        phi(iv.first.c[0] + 1, iv.first.c[1] + j, iv.first.c[2] + 0) -
        phi(iv.first.c[0] + 0, iv.first.c[1] + j, iv.first.c[2] + 1) +
        phi(iv.first.c[0] + 0, iv.first.c[1] + j, iv.first.c[2] + 0);
    HA(2, 0) += lambdaA[j][1] * h;
    HA(0, 2) += lambdaA[j][1] * h;
    HB(2, 0) += lambdaB[j][1] * h;
    HB(0, 2) += lambdaB[j][1] * h;
  }
  for (unsigned k = 0; k < 2; ++k) {
    const Scalar h =
        phi(iv.first.c[0] + 1, iv.first.c[1] + 1, iv.first.c[2] + k) -
        phi(iv.first.c[0] + 0, iv.first.c[1] + 1, iv.first.c[2] + k) -
        phi(iv.first.c[0] + 1, iv.first.c[1] + 0, iv.first.c[2] + k) +
        phi(iv.first.c[0] + 0, iv.first.c[1] + 0, iv.first.c[2] + k);
    HA(0, 1) += lambdaA[k][2] * h;
    HA(1, 0) += lambdaA[k][2] * h;
    HB(0, 1) += lambdaB[k][2] * h;
    HB(1, 0) += lambdaB[k][2] * h;
  }
}

// Explicit template instantiation
template class InterpolatedLevelSet<LEVELSET_INTERPOLATION_DEGREE>;

}  // namespace strandsim
