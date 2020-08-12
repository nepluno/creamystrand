/**
 * \copyright 2013 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#include "ProximityCollision.hh"

#include "../Core/ElasticStrand.hh"
#include "../Core/ElasticStrandUtils.hh"
#include "../Render/OpenGLHeaders.hh"
#include "../Utils/MathUtilities.hh"
#include "../Utils/TextLog.hh"

namespace strandsim {

const unsigned ProximityCollisionDatabase::s_maxAge = 3;

void ProximityCollision::generateTransformationMatrix() {
  Mat3x& E = transformationMatrix;
  E.col(0) = normal;

  Vec3x t1(-normal[1], normal[0], 0);
  const Scalar nt1 = t1.squaredNorm();

  if (isSmall(nt1)) {
    t1 = Vec3x(0, -normal[2], normal[1]);
    t1.normalize();
  } else {
    t1 /= std::sqrt(nt1);
  }

  E.col(1) = t1;
  E.col(2) = normal.cross(t1).normalized();
}

void ProximityCollision::updateTransformationMatrix(const Mat3x& previous) {
  Mat3x& E = transformationMatrix;
  E.col(0) = normal;

  Vec3x t1 =
      orthonormalParallelTransport(previous.col(1), previous.col(0), normal);
  orthoNormalize(t1, normal);
  E.col(1) = t1;

  E.col(2) = normal.cross(E.col(1)).normalized();
}

const ProximityCollision* ProximityCollisionDatabase::find(
    const ProximityCollision& needle) const {
  const Table& table = m_base[needle.objects.first.globalIndex];
  auto tableIt = table.find(needle.objects.second.globalIndex);

  const RecordKey key(needle.objects.first.vertex,
                      needle.objects.second.vertex);

#pragma omp atomic
  ++m_nQueries;

  if (tableIt != table.end()) {
    const Records& records = tableIt->second;

    auto recordIt = records.find(key);
    if (recordIt != records.end() && !recordIt->second.firstTime) {
#pragma omp atomic
      ++m_nFound;
      return &recordIt->second.collision;
    }
  }
  return NULL;
}

unsigned long ProximityCollision::computeSizeInBytes() const {
  unsigned long total = 0;

  // ContinuousTimeCollision* m_originalCTCollision;
  // std::cout << m_originalCTCollision << std::endl;
  total += sizeof(ContinuousTimeCollision*);
  // Vec3x normal;
  total += normal.size() * sizeof(Vec3x::Scalar);
  // Vec3x force;
  total += force.size() * sizeof(Vec3x::Scalar);
  // Scalar mu;
  total += sizeof(Scalar);
  // Scalar adhesion;
  total += sizeof(Scalar);
  // Scalar distance ;
  total += sizeof(Scalar);
  // Mat3x transformationMatrix ;
  total += transformationMatrix.size() * sizeof(Mat3x::Scalar);
  // std::pair<Object, Object> objects;
  total +=
      objects.first.computeSizeInBytes() + objects.second.computeSizeInBytes();

  return total;
}

void ProximityCollisionDatabase::insert(const ProximityCollision& collision) {
  Table& table = m_base[collision.objects.first.globalIndex];
  Records& records = table[collision.objects.second.globalIndex];

  const RecordKey key(collision.objects.first.vertex,
                      collision.objects.second.vertex);

  std::pair<Records::iterator, bool> res =
      records.insert(std::make_pair(key, Record()));
  Record& record = res.first->second;
  record.firstTime = res.second;
  record.age = s_maxAge;
  record.collision = collision;
}

int ProximityCollisionDatabase::numCollisions(int sIdx, int eIdx) const {
  auto itr = m_counter.find(RecordKey(sIdx, eIdx));
  if (itr != m_counter.end()) {
    return itr->second;
  }

  return 0;
}

bool ProximityCollisionDatabase::connectedLastStep(int sP, int iP, int sQ,
                                                   int iQ) const {
  if (sP == -1 || (sQ != -1 && sQ < sP)) {
    std::swap(sP, sQ);
    std::swap(iP, iQ);
  }

  const Table& table = m_base[sP];
  auto tableIt = table.find(sQ);
  if (tableIt == table.end()) return false;

  const RecordKey key(iP, iQ);
  const Records& records = tableIt->second;

  auto recordIt = records.find(key);
  return (recordIt != records.end() && recordIt->second.age == s_maxAge - 1);
}

void ProximityCollisionDatabase::ageAll() {
  unsigned numErased = 0;
  unsigned numRemaining = 0;

  // clear counter
  for (auto& counterPair : m_counter) {
    counterPair.second = 0;
  }

  for (unsigned i = 0; i < m_base.size(); ++i) {
    for (auto tableIt = m_base[i].begin(); tableIt != m_base[i].end();
         ++tableIt) {
      Records& records = tableIt->second;
      for (auto recordIt = records.begin(); recordIt != records.end();
           ++recordIt) {
        if (recordIt->second.age == s_maxAge) {
          ++m_counter[RecordKey(i, recordIt->first.first)];
          ++m_counter[RecordKey(tableIt->first, recordIt->first.second)];
        }
      }
    }
  }

  for (unsigned i = 0; i < m_base.size(); ++i) {
    for (auto tableIt = m_base[i].begin(); tableIt != m_base[i].end();
         ++tableIt) {
      Records& records = tableIt->second;
      for (auto recordIt = records.begin(); recordIt != records.end();) {
        if (--(recordIt->second.age)) {
          ++numRemaining;
          ++recordIt;
        } else {
          records.erase(recordIt++);
          ++numErased;
        }
      }
    }
  }

  DebugStream(g_log, "") << " Col database: " << numErased << " erased, "
                         << numRemaining << " remaining ";
  if (m_nQueries) {
    DebugStream(g_log, "") << " At last step: " << m_nQueries << " queries, "
                           << m_nFound << " found ( "
                           << ((100.f * m_nFound) / m_nQueries) << " % ) ";
  }

  m_nQueries = m_nFound = 0;
}

void ProximityCollisionDatabase::draw(
    const std::vector<ElasticStrand*>& strands) const {
  std::vector<float> normals;
  normals.reserve(6 * m_nFound);

  std::vector<float> colors;
  normals.reserve(6 * m_nFound);

  std::vector<float> forces;
  forces.reserve(6 * m_nFound);

  const Vec3f self =
      Vec3f(125, 195, 255) / 255.f;  // Hair/hair proximity, light blue
  const Vec3f prox =
      Vec3f(155, 250, 100) / 255.f;  // Hair/mesh proximity, brightgreen
  const Vec3f ct =
      Vec3f(100, 250, 175) / 255.f;  // Hair/mesh cont. time, turquoise
  const Vec3f inac = Vec3f(125, 125, 125) / 255.f;  // Inactive, middle grey
  const Vec3f force = Vec3f(255, 0, 255) / 255.f;   // Force, purple

  for (unsigned i = 0; i < m_base.size(); ++i) {
    for (auto tableIt = m_base[i].begin(); tableIt != m_base[i].end();
         ++tableIt) {
      const Records& records = tableIt->second;
      for (auto recordIt = records.begin(); recordIt != records.end();
           ++recordIt) {
        const bool active = (recordIt->second.age == s_maxAge);
        const ProximityCollision& c = recordIt->second.collision;

        Vec3x P = strands[c.objects.first.globalIndex]->getVertex(
            c.objects.first.vertex);
        if (c.objects.first.abscissa > 0) {
          P += c.objects.first.abscissa *
               (strands[c.objects.first.globalIndex]->getVertex(
                    c.objects.first.vertex + 1) -
                P);
        }

        normals.push_back(P[0]);
        normals.push_back(P[1]);
        normals.push_back(P[2]);

        if (active) {
          const Scalar n2Force = c.force.squaredNorm();
          if (!isSmall(n2Force)) {
            const Scalar len =
                1. - std::log(n2Force) / std::log(SMALL_NUMBER<Scalar>());
            const Vec3x f = c.force * len / std::sqrt(n2Force);

            forces.push_back(P[0]);
            forces.push_back(P[1]);
            forces.push_back(P[2]);

            const Vec3x Pf = P + c.transformationMatrix * f;

            forces.push_back(Pf[0]);
            forces.push_back(Pf[1]);
            forces.push_back(Pf[2]);
          }
        }

        P += c.normal;

        normals.push_back(P[0]);
        normals.push_back(P[1]);
        normals.push_back(P[2]);

        const Vec3f& color = active
                                 ? ((c.objects.second.globalIndex == -1)
                                        ? (c.m_originalCTCollision ? ct : prox)
                                        : self)
                                 : inac;

        colors.push_back(color[0]);
        colors.push_back(color[1]);
        colors.push_back(color[2]);
        colors.push_back(color[0]);
        colors.push_back(color[1]);
        colors.push_back(color[2]);
      }
    }
  }

  glPointSize(3.f);
  glLineWidth(1.f);

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);

  glVertexPointer(3, GL_FLOAT, 0, &normals[0]);
  glColorPointer(3, GL_FLOAT, 0, &colors[0]);

  glDrawArrays(GL_LINES, 0, normals.size() / 3);
  //    glDrawArrays( GL_POINTS, 0, normals.size() / 3 );

  glDisableClientState(GL_COLOR_ARRAY);

  glColor3f(force[0], force[1], force[2]);
  glVertexPointer(3, GL_FLOAT, 0, &forces[0]);
  glDrawArrays(GL_LINES, 0, forces.size() / 3);

  glDisableClientState(GL_VERTEX_ARRAY);
}

unsigned long ProximityCollisionDatabase::computeSizeInBytes() const {
  unsigned long total = 0;

  // static const unsigned s_maxAge ;
  total += sizeof(unsigned);
  // mutable unsigned m_nQueries;
  total += sizeof(unsigned);
  //  mutable unsigned m_nFound;
  total += sizeof(unsigned);
  // Base m_base : std::vector< std::unordered_map<int, std::map< std::pair<
  // int, int >, Record > > > Base;
  for (Base::size_type i = 0; i < m_base.size(); ++i) {
    // Each entry of m_base is a hash map: std::unordered_map<int, Records>
    total += sizeof(int);
    // For each entry of the hash map
    for (Table::const_iterator itr = m_base[i].begin(); itr != m_base[i].end();
         ++itr) {
      // Each Records is a map of type std::map< RecordKey, Record > Records
      // RecordKey is two ints
      total += 2 * sizeof(int);
      for (Records::const_iterator itr2 = (*itr).second.begin();
           itr2 != (*itr).second.end(); ++itr2) {
        // struct Record
        //{
        //  ProximityCollision collision ;
        //  unsigned short age ;
        //  bool firstTime ;
        //} ;
        total += sizeof(unsigned short);
        total += sizeof(bool);
        total += (*itr2).second.collision.computeSizeInBytes();
      }
    }
  }

  for (auto itr = m_counter.begin(); itr != m_counter.end(); ++itr) {
    // Each Records is a map of type std::map< RecordKey, int > Records
    // RecordKey is two ints, the counter is one int
    total += 3 * sizeof(int);
  }

  return total;
}

}  // namespace strandsim
