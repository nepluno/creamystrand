/**
 * \copyright 2012 Gilles Daviet
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef S_LEVELSET_H_
#define S_LEVELSET_H_

#include <fstream>
#include <string>
#include <vector>

#include "../Collision/TriangularMesh.hh"
#include "../Core/Definitions.hh"
#include "../Utils/ThreadUtils.hh"
#include "Bridson/array3.hh"
#include "Bridson/util.hh"
#include "LevelSetFwd.hh"

namespace strandsim {

namespace bridson {
typedef Array3<Scalar, Array1<Scalar> > Array3x;
}

typedef std::vector<unsigned int> Indices;
typedef std::vector<Vec3i> Vec3Indices;

template <int deg>
class InterpolatedLevelSet {
 public:
  struct AsyncDataPack {
    std::string m_szfn;
    Vec3x m_origin;
    Scalar m_dx;
    bridson::Array3x m_phi;
    bridson::Array3<Vec3x, bridson::Array1<Vec3x> > m_phiVel;
  };

  struct Coord {
    // Ok, I should probably rename those
    Vec3x f;            // point (in LS frame)
    Eigen::Vector3i c;  // cell
    Vec3x l;            // interp coefficients inside cell
    Scalar t;           // interp coeff on an edge
    Scalar d;           // LS value at this point
  };
  typedef std::pair<Coord, Coord> Interval;

  InterpolatedLevelSet();

  ~InterpolatedLevelSet();

  static void calculateLevelSetSize(const FaceArray &triangles,
                                    const Indices &triIndices,
                                    const Vec3xArray &x, const Vec3xArray &v,
                                    Vec3x &origin, Vec3x &center);

  void buildLevelSet(const FaceArray &triangles, const Indices &triIndices,
                     const Vec3xArray &x, const Vec3xArray &v,
                     const Vec3x &origin, const Scalar dx, const int nx,
                     const int ny, const int nz,
                     const Mat4x &transformMatrix = Mat4x::Identity());

  void buildLevelSetFromFile(std::ifstream &levelSetFile,
                             const FaceArray &triangles,
                             const Indices &triIndices, const Vec3xArray &x,
                             const Vec3xArray &v, const Vec3x &origin,
                             const Scalar dx, const int nx, const int ny,
                             const int nz,
                             const Mat4x &transformMatrix = Mat4x::Identity());

  void buildAdaptativeLevelSet(
      const FaceArray &triangles, const Indices &triIndices,
      const Vec3xArray &x, const Vec3xArray &v,
      const Mat4x &transformMatrix = Mat4x::Identity());

  void update();

  void updateFromFile(std::ifstream &levelSetFile);

  Scalar getLevelSetValue(const Vec3x &x) const;
  Scalar getLevelSetValue(int i, int j, int k) const;
  Scalar getEdgeValue(const Vec3x &start, const Vec3x &end, Coord &min) const;
  Scalar getIntegralValue(const Interval &iv) const;
  Vec3x getWorldSpacePos(int i, int j, int k) const;

  void getInteriorPoints(const Vec3x &start, const Vec3x &end,
                         const Scalar thickness, std::vector<Coord> &points,
                         bool stopAtFirstMin = false) const;
  void getInteriorIntervals(const Vec3x &start, const Vec3x &end,
                            const Scalar thickness,
                            std::vector<Interval> &intervals) const;

  Scalar getLevelSetValueVelocity(const Vec3x &x, Vec3x &v) const;

  void getGradient(const Vec3x &x, Vec3x &grad, bool rotated = true) const;
  void getEdgeGradient(const Vec3x &start, const Vec3x &end, const Coord &min,
                       Vec3x &grad, int *gradLambdaAxis = NULL) const;
  void getIntegralGradient(const Interval &iv, Vec3x &gradStart,
                           Vec3x &gradEnd) const;

  void getHessian(const Vec3x &x, Mat3x &hess, bool rotated = true) const;
  void getEdgeHessian(const Vec3x &start, const Vec3x &end, const Scalar lambda,
                      const int gradLambdaAxis,
                      Eigen::Matrix<Scalar, 6, 6> &hess) const;
  void getIntegralHessian(const Interval &iv, Mat6x &hess) const;

  unsigned getClosestTriangle(const Vec3x &x, Scalar &t1, Scalar &t2,
                              Scalar &t3) const;
  Vec3x getCurrentPosition(unsigned t, Scalar t1, Scalar t2, Scalar t3) const;

  const Vec3x &getOrigin() const { return m_origin; }

  Scalar getGridSize() const { return m_dx; }

  int getNbrX() const { return m_phi.ni; }

  int getNbrY() const { return m_phi.nj; }

  int getNbrZ() const { return m_phi.nk; }

  Vec3x getDims() const {
    return (m_length + 2 * (m_origin - m_unpaddedOrigin));
  }

  const bridson::Array3x &getPhi() const { return m_phi; }

  const bridson::Array3<Vec3x, bridson::Array1<Vec3x> > &getPhiVel() const {
    return m_phiVel;
  }

  void draw(const Scalar scale = 1.0) const;
  void drawInterior(const Scalar scale = 1.0, Scalar centerValueToDraw = 0.0f,
                    Scalar thicknessToDraw = 0.01f);

  void writeFileAsync(const std::string &szfn) const;

  void writeFile(std::ofstream &levelSetFile) const;

  void loadFile(std::ifstream &levelSetFile);

  size_t getFileSize(const Vec3x &origin, const Vec3x &center,
                     const Scalar &dx) const;

  bool isInitialized() const { return m_initialized; }

  // This stores the transformation matrix for the mesh the level set is created
  // from
  void setTransformationMatrix(const Mat4x &i_matrix);
  const Mat4x &getTransformationMatrix() const {
    return m_relativeTransformMatrix;
  }
  const Mat4x &getInverseTransformationMatrix() const {
    return m_invRelativeTransformMatrix;
  }

  Mat4x getOriginTransformationMatrix() const {
    Mat4x otm = getTransformationMatrix();
    otm.block<3, 1>(0, 3) -= m_unpaddedOrigin;
    return otm;
  }

  void setScale(const Scalar scale);

  bool intersectsAABB(const std::pair<Vec3x, Vec3x> &aaBB) const;

  // For debug purposes only
  void set(bridson::Array3x phi, Vec3x origin, Vec3x center, Scalar dx);

  void computeAABB();

  void computeBarycentric(Coord &c) const;
  void normalizeInterval(Interval &i) const;

  Scalar getInterpolated(int i, int j, int k, Scalar fi, Scalar fj,
                         Scalar fk) const;
  inline Scalar getInterpolated(const Coord &c) const;
  inline Scalar getTrilerp(const Coord &c) const;

  Scalar getFractionInsideU(const Vec3x &x) const;
  Scalar getFractionInsideV(const Vec3x &x) const;
  Scalar getFractionInsideW(const Vec3x &x) const;
  Scalar getFractionInsideU(int i, int j, int k, Scalar fi, Scalar fj,
                            Scalar fk) const;
  Scalar getFractionInsideV(int i, int j, int k, Scalar fi, Scalar fj,
                            Scalar fk) const;
  Scalar getFractionInsideW(int i, int j, int k, Scalar fi, Scalar fj,
                            Scalar fk) const;
  Scalar FractionInside(Scalar phi_bl, Scalar phi_br, Scalar phi_tl,
                        Scalar phi_tr) const;
  Scalar FractionInside(Scalar phi_left, Scalar phi_right) const;
  Vec3x getInterpolatedGradient(int i, int j, int k, Scalar fi, Scalar fj,
                                Scalar fk) const;
  Mat3x getInterpolatedHessian(int i, int j, int k, Scalar fi, Scalar fj,
                               Scalar fk) const;

  void getEdgeSamplingPoints(const Vec3x &start, const Vec3x &end, Coord &min,
                             std::vector<Coord> *points = NULL,
                             Scalar maxVal = 0.,
                             bool stopAtFirstMin = false) const;
  void getEdgeIntervals(const Vec3x &start, const Vec3x &end,
                        const Scalar maxVal,
                        std::vector<Interval> &intervals) const;

 protected:
  bridson::Array3x m_phi;
  bridson::Array3<Vec3x, bridson::Array1<Vec3x> > m_phiVel;
  Vec3x m_origin;
  Vec3x m_length;
  Scalar m_dx;
  Scalar m_scale;
  Vec3x m_aaBBMin;
  Vec3x m_aaBBMax;

  bool m_initialized;
  volatile bool m_needsRebuilding;
  bool m_adaptativeGridSize;
  volatile bool m_transformChanged;

  mutable MutexType m_levelSetMutex;
  mutable MutexType m_transformMutex;

  // Store information about the mesh, in order to computes points displacements
  bridson::Array3i m_closest_tri;

  // Store parameters for asynchronous creation of level set
  Vec3x m_newOrigin;
  Vec3x m_newCenter;

  Vec3x m_unpaddedOrigin;

  const Vec3xArray *m_x;
  const Vec3xArray *m_v;
  const Indices *m_triangles;
  const FaceArray *m_triIndices;

  // The mesh may not be at the origin when it is created. So remember where it
  // was so that later when we're asked to sample the level set we can remove
  // the initial transform before working out the true place to ask.
  Mat4x m_transformMatrixAtCreation;
  Mat4x m_currentTransformMatrix;
  Mat4x m_relativeTransformMatrix;
  Mat4x m_invRelativeTransformMatrix;

  // Debug data to check whether the values returned work for meshes that move
  // or were created with a transform on them
  std::vector<Vec4x> m_realRequestPositions;
  std::vector<Vec4x> m_transformedRequestPositions;
  std::vector<Vec4x> m_grad;
  std::vector<Vec4x> m_gradPosition;
};

/*
 * Utility functions
 */

void computeLevelSet(const FaceArray &triangles, const Indices &triIndices,
                     const Vec3xArray &x, const Vec3xArray &v, Vec3x &origin,
                     const Vec3x &center, Scalar &dx, bridson::Array3x &phi,
                     bridson::Array3<Vec3x, bridson::Array1<Vec3x> > &phiVel,
                     bridson::Array3i &closest_tri, bool adaptativeDx);

void sweep(const FaceArray &tri, const Vec3xArray &x, const Vec3xArray &v,
           bridson::Array3x &phi,
           bridson::Array3<Vec3x, bridson::Array1<Vec3x> > &phi_vel,
           bridson::Array3i &closest_tri, const Vec3x &origin, const Scalar dx,
           const int di, const int dj, const int dk);

int orientation(Scalar x1, Scalar y1, Scalar x2, Scalar y2,
                Scalar &twice_signed_area);

bool point_in_triangle_2d(const Scalar x0, const Scalar y0, Scalar x1,
                          Scalar y1, Scalar x2, Scalar y2, Scalar x3, Scalar y3,
                          Scalar &a, Scalar &b, Scalar &c);

Scalar point_triangle_distance(const Vec3x &p, const Vec3x &a, const Vec3x &b,
                               const Vec3x &c, Scalar &t1, Scalar &t2,
                               Scalar &t3);

void check_neighbour(const FaceArray &tri, const Vec3xArray &x,
                     const Vec3xArray &v, bridson::Array3x &phi,
                     bridson::Array3<Vec3x, bridson::Array1<Vec3x> > &phi_vel,
                     bridson::Array3i &closest_tri, const Vec3x &gx, int i0,
                     int j0, int k0, int i1, int j1, int k1);

static void computeTrilinearValueCoeffs(
    const bridson::Array3x &phi, const InterpolatedLevelSet<1>::Interval &iv,
    Vec4x &coeffs);
static void computeBilinearGradientCoeffs(
    const bridson::Array3x &phi, const InterpolatedLevelSet<1>::Interval &iv,
    Mat3x &coeffs);
static void computeLinearHessianCoeffs(
    const bridson::Array3x &phi, const InterpolatedLevelSet<1>::Interval &iv,
    Mat3x &HA, Mat3x &HB);

}  // namespace strandsim

#endif
