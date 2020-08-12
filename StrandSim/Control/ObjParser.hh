/**
 * \copyright 2009 Breannan Smith
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef OBJPARSER_HH
#define OBJPARSER_HH

#include <string>

#include "../Collision/TriangularMesh.hh"
#include "../Forces/Bridson/array3.hh"

namespace strandsim {
/**
 * Loades meshes from OBJ format files.
 */
class ObjParser {
 public:
  ObjParser(){};
  ObjParser(const std::string& obj_file_name, TriangularMesh* mesh,
            bool inverted);
  /**
   * Loads a triangle mesh into tri_mesh from OBJ format file obj_file_name.
   *
   * \param[in] obj_file_name String filename of OBJ format file.
   * \param[in] tri_mesh Empty tirangle mesh to load OBJ file into.
   */
  bool loadTriangularMesh(const std::string& obj_file_name,
                          TriangularMesh& tri_mesh, bool inverted = false);

 private:
  void parsevCommand(std::istringstream& commandstream,
                     TriangularMesh& tri_mesh);
  void parsetrianglefCommand(std::istringstream& commandstream,
                             TriangularMesh& tri_mesh);
  void tokenize(const std::string& str, std::vector<std::string>& tokens,
                const std::string& delimiters = " ");

  bool m_vn_vt_unsupported_printed;
  bool m_vt_unsupported_printed;
  bool m_vn_unsupported_printed;
  bool m_o_unsupported_printed;
  bool m_g_unsupported_printed;
  bool m_mtllib_unsupported;
  bool m_usemtl_unsupported;
  bool m_s_unsupported;

  bool m_successful_parse;

  int m_line_num;
};

/**
 * Loads points from COFF format files.
 */
class CoffParser {
 public:
  bool loadParticles(const std::string& coff_file_name,
                     std::vector<Vec3x>& particles);
};

template <typename S, typename T = S>
class ArrayWriter {
 public:
  bool writeArray3(const std::string& array_file_name,
                   const bridson::Array3<S, bridson::Array1<S> >& arr);
  bool writeArray3(const std::string& array_file_name,
                   const std::vector<S>& arr, int ni, int nj, int nk);
};
}  // namespace strandsim

#endif
