/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "ReferenceFrames.hh"

#include "../Core/ElasticStrandUtils.hh"
#include "../Utils/MathUtilities.hh"

namespace strandsim {

void ReferenceFrames1::storeInitialFrames(const Vec3x& initRefFrame1) {
  m_value.resize(m_size);
  const Vec3xArray& tangents = m_tangents.get();
  const Vec3x& tan0 = tangents[0];
  assert(isApproxUnit(tan0));

  // Do we have an even approximately valid initial reference frame?
  if (initRefFrame1.squaredNorm() > 0.5 &&
      fabs(initRefFrame1.dot(tan0)) < 0.25) {
    // If so, just project it on the plane normal to the tangent vector
    const Vec3x projectedInitRefFrame1 =
        (initRefFrame1 - initRefFrame1.dot(tan0) * tan0).normalized();
    m_value[0] = projectedInitRefFrame1;
  } else  // If a valid initial first reference frame hasn't been provided, use
          // an arbitrary one
  {
    m_value[0] = findNormal<3>(tan0);
  }

  // Next initial reference frames are obtained by space-parallel transportation
  // along the rod
  for (IndexType vtx = 1; vtx < size(); ++vtx) {
    m_value[vtx] = orthonormalParallelTransport(
        m_value[vtx - 1], tangents[vtx - 1], tangents[vtx]);
    orthoNormalize(m_value[vtx], tangents[vtx]);
  }

  // Store tangents backup for time-parallel transport
  m_previousTangents = tangents;

  setClean();
  setDependentsDirty();
}

void ReferenceFrames1::compute() {
  m_value.resize(m_size);
  const Vec3xArray& tangents = m_tangents.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    Vec3x& previousTangent = m_previousTangents[vtx];
    const Vec3x& currentTangent = tangents[vtx];

    m_value[vtx] = orthonormalParallelTransport(m_value[vtx], previousTangent,
                                                currentTangent);
    orthoNormalize(m_value[vtx], currentTangent);

    // Store the current tangent for the next time-parallel transportation
    previousTangent = currentTangent;
  }
  setDependentsDirty();
}

bool ReferenceFrames1::checkNormality() {
  bool normal = true;

  const Vec3xArray& tangents = m_previousTangents;

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    if (!isSmall(m_value[vtx].dot(tangents[vtx]))) {
      normal = false;
      ErrorStream(g_log, "")
          << "ReferenceFrames1::checkNormality() fails at vtx = " << vtx;
      break;
    }
  }
  return normal;
}

void ReferenceFrames2::compute() {
  m_value.resize(m_size);
  const Vec3xArray& tangents = m_tangents.get();
  const Vec3xArray& referenceFrames1 = m_referenceFrames1.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    m_value[vtx] = tangents[vtx].cross(referenceFrames1[vtx]);
  }

  setDependentsDirty();
}

void ReferenceTwists::compute() {
  m_value.resize(m_size);
  const Vec3xArray& tangents = m_tangents.get();
  const Vec3xArray& referenceFrames1 = m_referenceFrames1.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    const Vec3x& u0 = referenceFrames1[vtx - 1];
    const Vec3x& u1 = referenceFrames1[vtx];
    const Vec3x& tangent = tangents[vtx];

    // transport reference frame to next edge
    Vec3x ut = orthonormalParallelTransport(u0, tangents[vtx - 1], tangent);

    // rotate by current value of reference twist
    const Scalar beforeTwist = m_value[vtx];
    rotateAxisAngle(ut, tangent, beforeTwist);

    // compute increment to reference twist to align reference frames
    m_value[vtx] = beforeTwist + signedAngle(ut, u1, tangent);
  }

  setDependentsDirty();
}

}  // namespace strandsim
