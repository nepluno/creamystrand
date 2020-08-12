/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ELEMENTPROXY_HH_
#define ELEMENTPROXY_HH_

#include "../Core/BoundingBox.hh"
#include "../Core/ElasticStrand.hh"
#include "../Dynamic/MeshScriptingController.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "TriangularMesh.hh"

namespace strandsim {

typedef BoundingBox<float> BBoxType;

class ElementProxy {
 public:
  ElementProxy() {}
  virtual ~ElementProxy() {}

  virtual void computeBoundingBox(BBoxType& boundingBox,
                                  bool statique) const = 0;

  void updateBoundingBox(bool statique) {
    computeBoundingBox(m_boundingBox, statique);
  }

  const BBoxType& getBoundingBox() const { return m_boundingBox; }

  void resetBoundingBox() { m_boundingBox.reset(); }

  friend std::ostream& operator<<(std::ostream& os, const ElementProxy& elem);

 protected:
  virtual void print(std::ostream& os) const = 0;

  BBoxType m_boundingBox;
};

static const int FIRST_NON_IMMUNE = 2;

class EdgeProxy : public ElementProxy {
 public:
  EdgeProxy(ElasticStrand& strand, int vertexIndex)
      : m_strand(strand), m_vertexIndex(vertexIndex) {}

  void computeBoundingBox(BBoxType& boundingBox, bool statique) const {
    boundingBox.reset();

    if (shouldIgnore()) return;

    const Vec3f& vtx0 = m_strand.getVertex(m_vertexIndex).cast<float>();
    const Vec3f& vtx1 = m_strand.getVertex(m_vertexIndex + 1).cast<float>();

    boundingBox.insert(
        vtx0 -
        m_strand.dynamics().getDisplacement(m_vertexIndex).cast<float>());
    boundingBox.insert(
        vtx1 -
        m_strand.dynamics().getDisplacement(m_vertexIndex + 1).cast<float>());

    if (!statique)  // If non-static, the strand's current contain the post-step
                    // positions
    {
      boundingBox.insert(vtx0);
      boundingBox.insert(vtx1);
    }
  }

  bool shouldIgnore() const {
    return !m_strand.activelySimulated() ||
           m_strand.dynamics().isImmune(getVertexIndex());
  }

  const ElasticStrand* getStrandPointer() const { return &m_strand; }
  ElasticStrand* getStrandPointer() { return &m_strand; }
  int getVertexIndex() const { return m_vertexIndex; }

  ElasticStrand& getStrand() { return m_strand; }

 protected:
  virtual void print(std::ostream& os) const;

  ElasticStrand& m_strand;  // Not const
  const int m_vertexIndex;
};

class CylinderProxy : public EdgeProxy {
 public:
  CylinderProxy(ElasticStrand& strand, int vertexIndex, Scalar radius)
      : EdgeProxy(strand, vertexIndex), m_radius(radius) {}

  void computeBoundingBox(BBoxType& boundingBox, bool statique) const {
    if (shouldIgnore()) {
      boundingBox.reset();
      return;
    }

    Vec3x min, max;
    m_strand.getFutureAABB(m_vertexIndex, min, max,
                           CollisionParameters::EXTERNAL);
    boundingBox.min = min.cast<float>();
    boundingBox.max = max.cast<float>();

    if (!statique)  // If non-static, the strand's current contain the post-step
                    // positions
    {
      m_strand.getAABB(m_vertexIndex, min, max, CollisionParameters::EXTERNAL);
      boundingBox.min = boundingBox.min.cwiseMin(min.cast<float>());
      boundingBox.max = boundingBox.max.cwiseMax(max.cast<float>());
    }
  }

 protected:
  const Scalar m_radius;
};

class FaceProxy : public ElementProxy {
 public:
  FaceProxy(const unsigned faceIndex,
            const std::shared_ptr<MeshScriptingController>& controller)
      : m_faceIndex(faceIndex),
        m_face(controller->getCurrentMesh()->getFace(faceIndex)),
        m_controller(controller) {}

  void computeBoundingBox(BBoxType& boundingBox, bool statique) const {
    boundingBox.reset();

    const Vec3f vtx0 = getVertex(0).cast<float>();
    const Vec3f vtx1 = getVertex(1).cast<float>();
    const Vec3f vtx2 = getVertex(2).cast<float>();

    boundingBox.insert(vtx0 - getDisplacement(0).cast<float>());
    boundingBox.insert(vtx1 - getDisplacement(1).cast<float>());
    boundingBox.insert(vtx2 - getDisplacement(2).cast<float>());

    if (!statique)  // If non-static, the strand's current contain the post-step
                    // positions
    {
      boundingBox.insert(vtx0);
      boundingBox.insert(vtx1);
      boundingBox.insert(vtx2);
    }
  }

  Vec3x getVertex(short apex) const {
    return m_controller->getCurrentMesh()->getVertex(m_face.idx[apex]);
  }

  Vec3x getDisplacement(short apex) const {
    return m_controller->getCurrentMesh()->getDisplacement(m_face.idx[apex]);
  }

  int getVertexIdx(short apex) const { return m_face.idx[apex]; }

  bool hasEnabledApex(short apex) const {
    const std::vector<bool>& enabledVertices =
        m_controller->getEnabledVertices();

    if (enabledVertices.empty()) {
      return true;
    } else {
      return enabledVertices[m_face.idx[apex]];
    }
  }

  bool allApicesEnabled() const {
    const std::vector<bool>& enabledVertices =
        m_controller->getEnabledVertices();

    if (enabledVertices.empty()) {
      return true;
    } else {
      return enabledVertices[m_face.idx[0]] && enabledVertices[m_face.idx[1]] &&
             enabledVertices[m_face.idx[2]];
    }
  }

  double getFrictionCoefficient(double b0, double b1, double b2) const {
    const std::vector<double>& frictionCoefficients =
        m_controller->getFrictionCoefficients();

    if (frictionCoefficients.empty()) {
      return m_controller->getDefaultFrictionCoefficient();
    } else {
      return b0 * frictionCoefficients[m_face.idx[0]] +
             b1 * frictionCoefficients[m_face.idx[1]] +
             b2 * frictionCoefficients[m_face.idx[2]];
    }
  }

  const std::shared_ptr<TriangularMesh> getMesh() const {
    return m_controller->getCurrentMesh();
  }

  Vec3x getNormal() const {
    const Vec3x& q0 = getVertex(0);
    const Vec3x& q1 = getVertex(1);
    const Vec3x& q2 = getVertex(2);

    return (q1 - q0).cross(q2 - q0).normalized();
  }

  bool collideOnBothSides() const { return m_controller->collideOnBothSides(); }

  short knowsNormalSign(bool atPreviousStep, const ElasticStrand& strand,
                        unsigned vertex) const {
    return m_controller->knowsNormalSign(atPreviousStep, m_faceIndex,
                                         strand.getGlobalIndex(), vertex);
  }

  void setNormalSign(short sign, float offset, const ElasticStrand& strand,
                     unsigned vertex) const {
    m_controller->setNormalSign(sign, offset, m_faceIndex,
                                strand.getGlobalIndex(), vertex);
  }

  const TriangularFace& getFace() const { return m_face; }

  //! Guarranteed to start with a zero and end with two zeros
  uint32_t uniqueId() const {
    // Not the real Id, but will do as long as it is unique
    // &m_face is aligned on 4 bytes, so
    // we can afford to lose the 4 least significant bits of the adress
    // Since we will zero-out the last two, we can still shift by two bits

    return ((((size_t)&m_face) >> 2) & 0x7FFFFFFCul);
  }

 protected:
  virtual void print(std::ostream& os) const;

  //    const TriangularMesh* m_mesh;
  //    const std::vector<bool>* m_enabledVertices;
  const unsigned m_faceIndex;
  const TriangularFace& m_face;
  std::shared_ptr<MeshScriptingController> m_controller;
};

class ElementProxyBBoxFunctor {
 public:
  ElementProxyBBoxFunctor(std::vector<ElementProxy*>& elements)
      : m_elements(elements) {
    for (unsigned i = 0; i < m_elements.size(); ++i) {
      if (m_elements[i]->getBoundingBox().isValid()) {
        m_valid.push_back(i);
      }
    }
    for (unsigned i = 0; i < m_valid.size(); ++i) {
      swap(i, m_valid[i]);
    }
  }

  unsigned int size() const { return (unsigned int)m_valid.size(); }

  const BBoxType& operator[](const unsigned int i) const {
    return m_elements[i]->getBoundingBox();
  }

  void swap(unsigned int i, unsigned int j) {
    std::swap(m_elements[i], m_elements[j]);
  }

 private:
  std::vector<ElementProxy*>& m_elements;
  std::vector<unsigned> m_valid;
};

class ElementProxySortedAABBFunctor {
 public:
  ElementProxySortedAABBFunctor(const std::vector<ElementProxy*>& proxies) {
    for (unsigned i = 0; i < proxies.size(); ++i) {
      if (dynamic_cast<EdgeProxy*>(proxies[i])) {
        m_edgeProxies.push_back(proxies[i]);
      } else {
        m_otherProxies.push_back(proxies[i]);
      }
    }
  }

  ElementProxy* elementProxy(unsigned elementID) const {
    const unsigned n = numEdgeProxies();
    return elementID < n ? m_edgeProxies[elementID]
                         : m_otherProxies[elementID - n];
  }

  void getAABB(unsigned elementID, Vec3x& min, Vec3x& max) const {
    const BBoxType& bbox = elementProxy(elementID)->getBoundingBox();
    min = bbox.min.cast<Scalar>();
    max = bbox.max.cast<Scalar>();
  }

  unsigned numEdgeProxies() const { return m_edgeProxies.size(); }

 private:
  std::vector<ElementProxy*> m_edgeProxies;
  std::vector<ElementProxy*> m_otherProxies;
};

} /* namespace strandsim */
#endif /* ELEMENTPROXY_HH_ */
