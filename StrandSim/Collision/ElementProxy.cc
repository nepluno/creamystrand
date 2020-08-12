/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "ElementProxy.hh"

#include "../Dynamic/MeshScriptingController.hh"
#include "BVH.hh"

namespace strandsim {

std::ostream& operator<<(std::ostream& os, const ElementProxy& elem) {
  elem.print(os);

  return os;
}

void EdgeProxy::print(std::ostream& os) const {
  os << "edge: " << &m_strand << ' ' << m_vertexIndex << ": "
     << m_strand.getVertex(m_vertexIndex).format(EIGEN_VECTOR_IO_FORMAT)
     << " --- "
     << m_strand.getVertex(m_vertexIndex + 1).format(EIGEN_VECTOR_IO_FORMAT);
}

void FaceProxy::print(std::ostream& os) const {
  auto mesh = getMesh();

  os << "triangle: " << mesh << '\n';

  os << "past position:\n";
  os << m_face.idx[0] << ": "
     << (mesh->getVertex(m_face.idx[0]) - mesh->getDisplacement(m_face.idx[0]))
            .format(EIGEN_VECTOR_IO_FORMAT)
     << '\n';
  os << m_face.idx[1] << ": "
     << (mesh->getVertex(m_face.idx[1]) - mesh->getDisplacement(m_face.idx[1]))
            .format(EIGEN_VECTOR_IO_FORMAT)
     << '\n';
  os << m_face.idx[2] << ": "
     << (mesh->getVertex(m_face.idx[2]) - mesh->getDisplacement(m_face.idx[2]))
            .format(EIGEN_VECTOR_IO_FORMAT)
     << '\n';

  os << "current position:\n";
  os << m_face.idx[0] << ": "
     << mesh->getVertex(m_face.idx[0]).format(EIGEN_VECTOR_IO_FORMAT) << '\n';
  os << m_face.idx[1] << ": "
     << mesh->getVertex(m_face.idx[1]).format(EIGEN_VECTOR_IO_FORMAT) << '\n';
  os << m_face.idx[2] << ": "
     << mesh->getVertex(m_face.idx[2]).format(EIGEN_VECTOR_IO_FORMAT);
}

} /* namespace strandsim */
