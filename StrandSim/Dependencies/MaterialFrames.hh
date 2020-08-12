/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MATERIALFRAMES_HH_
#define MATERIALFRAMES_HH_

#include "ReferenceFrames.hh"

namespace strandsim {

/**
 * Unit: no dimension
 */
template <int FrameN>
class MaterialFrames : public DependencyNode<Vec3xArray> {
 public:
  MaterialFrames(TrigThetas& trigThetas, ReferenceFrames1& referenceFrames1,
                 ReferenceFrames2& referenceFrames2)
      : DependencyNode<Vec3xArray>(0, referenceFrames1.size()),
        m_trigThetas(trigThetas),
        m_referenceFrames1(referenceFrames1),
        m_referenceFrames2(referenceFrames2) {
    m_trigThetas.addDependent(this);
    m_referenceFrames1.addDependent(this);
    m_referenceFrames2.addDependent(this);
  }

  virtual const char* name() const;

 protected:
  virtual void compute();
  Vec3x linearMix(const Vec3x& u, const Vec3x& v, Scalar s, Scalar c);

  TrigThetas& m_trigThetas;
  ReferenceFrames1& m_referenceFrames1;
  ReferenceFrames2& m_referenceFrames2;
};

}  // namespace strandsim

#endif /* MATERIALFRAMES_HH_ */
