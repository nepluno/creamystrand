/**
 * \copyright 2009 Breannan Smith
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "ObjParser.hh"

#include <fstream>

namespace strandsim {
ObjParser::ObjParser(const std::string& obj_file_name, TriangularMesh* mesh,
                     bool inverted) {
  loadTriangularMesh(obj_file_name, *mesh, inverted);
}

bool ObjParser::loadTriangularMesh(const std::string& obj_file_name,
                                   TriangularMesh& tri_mesh, bool inverted) {
  m_vn_vt_unsupported_printed = false;
  m_vn_unsupported_printed = false;
  m_vt_unsupported_printed = false;
  m_o_unsupported_printed = false;
  m_g_unsupported_printed = false;
  m_mtllib_unsupported = false;
  m_usemtl_unsupported = false;
  m_s_unsupported = false;

  m_successful_parse = true;

  m_line_num = 0;

  std::ifstream obj_file;
  obj_file.open(obj_file_name.c_str());
  if (!obj_file.is_open()) {
    std::cout << "OBJPARSER: Failed to open file " << obj_file_name
              << std::endl;
    return false;
  }

  std::string obj_command;
  while (!obj_file.eof()) {
    ++m_line_num;

    // Read a single line at a time
    getline(obj_file, obj_command);

    // Use a string stream for easy tokenizing
    std::istringstream commandstream(obj_command);

    // First element of a command is the command's name
    std::string command;
    commandstream >> command;

    // Comment command
    if (command == "#") continue;
    // Vertex command
    else if (command == "v")
      parsevCommand(commandstream, tri_mesh);
    // Face command
    else if (command == "f")
      parsetrianglefCommand(commandstream, tri_mesh);
    // Unsupported commands
    else if (command == "vt") {
      // if( !m_vt_unsupported_printed ) std::cout << "OBJPARSER: vt command not
      // supported by obj parsesr." << std::endl;
      m_vt_unsupported_printed = true;
    } else if (command == "vn") {
      if (!m_vn_unsupported_printed)
        std::cout << "OBJPARSER: vn command not supported by obj parsesr."
                  << std::endl;
      m_vn_unsupported_printed = true;
    } else if (command == "o") {
      if (!m_o_unsupported_printed)
        std::cout << "OBJPARSER: o command not supported by obj parsesr."
                  << std::endl;
      m_o_unsupported_printed = true;
    } else if (command == "g") {
      if (!m_g_unsupported_printed)
        std::cout << "OBJPARSER: g command not supported by obj parsesr."
                  << std::endl;
      m_g_unsupported_printed = true;
    } else if (command == "mtllib") {
      if (!m_mtllib_unsupported)
        std::cout << "OBJPARSER: mtllib command not supported by obj parsesr."
                  << std::endl;
      m_mtllib_unsupported = true;
    } else if (command == "usemtl") {
      if (!m_usemtl_unsupported)
        std::cout << "OBJPARSER: usemtl command not supported by obj parsesr."
                  << std::endl;
      m_usemtl_unsupported = true;
    } else if (command == "s") {
      if (!m_s_unsupported)
        std::cout << "OBJPARSER: s command not supported by obj parsesr."
                  << std::endl;
      m_s_unsupported = true;
    }
  }
  tri_mesh.computeBBox();

  if (inverted) {
    for (TriangularFace& f : tri_mesh.m_faces) {
      f.invert();
    }
  }

  obj_file.close();
  return m_successful_parse;
}

void ObjParser::parsevCommand(std::istringstream& commandstream,
                              TriangularMesh& tri_mesh) {
  double x, y, z;
  if (!(commandstream >> x >> y >> z)) {
    std::cout << "OBJPARSER: Invalid vertex command on line " << m_line_num
              << ", vertex must have three real coordinates." << std::endl;
    m_successful_parse = false;
    return;
  }
  // Create a new vertex
  //    TriangularMesh::vertex_handle vh = tri_mesh.addVertex();
  //    tri_mesh.getVertex(vh) = Vec3d(x,y,z);
  tri_mesh.addVertex(Vec3d(x, y, z));
}

void ObjParser::parsetrianglefCommand(std::istringstream& commandstream,
                                      TriangularMesh& tri_mesh) {
  int xidx, yidx, zidx;

  std::vector<std::string> vertstrngs;
  std::string vertcmmnd;
  while (commandstream >> vertcmmnd) vertstrngs.push_back(vertcmmnd);

  if (vertstrngs.size() != 3) {
    std::cout << vertstrngs.size() << std::endl;
    std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
              << ", face must have three vertices." << std::endl;
    m_successful_parse = false;
    return;
  }

  std::string xstr = vertstrngs[0];
  std::string ystr = vertstrngs[1];
  std::string zstr = vertstrngs[2];

  // Attempt to extract vertex indices for x
  std::vector<std::string> vertex_spec;
  tokenize(xstr, vertex_spec, "/");
  if (vertex_spec.size() < 1 || vertex_spec.size() > 3) {
    std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
              << ", x vertex specification must be of form v[/vt][/vn]."
              << std::endl;
    m_successful_parse = false;
    return;
  }
  if (vertex_spec.size() != 1) {
    // if( !m_vn_vt_unsupported_printed ) std::cout << "OBJPARSER: Texture
    // coordinates and normals not supported by obj parsesr." << std::endl;
    m_vn_vt_unsupported_printed = true;
  }
  std::istringstream xstream(vertex_spec[0]);
  if (!(xstream >> xidx)) {
    std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
              << ", x vertex index must be an integer." << std::endl;
    m_successful_parse = false;
    return;
  }
  if (xidx > (int)tri_mesh.nv()) {
    std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
              << ", x vertex index greater than number of vertices."
              << std::endl;
    m_successful_parse = false;
    return;
  }
  if (xidx < 0)
    if (xidx + ((int)tri_mesh.nv()) < 0) {
      std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
                << ", relative x vertex index out of bounds." << std::endl;
      m_successful_parse = false;
      return;
    }
  // Attempt to extract vertex indices for y
  vertex_spec.clear();
  tokenize(ystr, vertex_spec, "/");
  if (vertex_spec.size() < 1 || vertex_spec.size() > 3) {
    std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
              << ", y vertex specification must be of form v[/vt][/vn]."
              << std::endl;
    m_successful_parse = false;
    return;
  }
  if (vertex_spec.size() != 1) {
    if (!m_vn_vt_unsupported_printed)
      std::cout << "OBJPARSER: Texture coordinates and normals not supported "
                   "by obj parsesr."
                << std::endl;
    m_vn_vt_unsupported_printed = true;
  }
  std::istringstream ystream(vertex_spec[0]);
  if (!(ystream >> yidx)) {
    std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
              << ", y vertex index must be an integer." << std::endl;
    m_successful_parse = false;
    return;
  }
  if (yidx > (int)tri_mesh.nv()) {
    std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
              << ", y vertex index greater than number of vertices."
              << std::endl;
    m_successful_parse = false;
    return;
  }
  if (yidx < 0)
    if (yidx + ((int)tri_mesh.nv()) < 0) {
      std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
                << ", relative y vertex index out of bounds." << std::endl;
      m_successful_parse = false;
      return;
    }
  // Attempt to extract vertex indices for z
  vertex_spec.clear();
  tokenize(zstr, vertex_spec, "/");
  if (vertex_spec.size() < 1 || vertex_spec.size() > 3) {
    std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
              << ", z vertex specification must be of form v[/vt][/vn]."
              << std::endl;
    m_successful_parse = false;
    return;
  }
  if (vertex_spec.size() != 1) {
    if (!m_vn_vt_unsupported_printed)
      std::cout << "OBJPARSER: Texture coordinates and normals not supported "
                   "by obj parsesr."
                << std::endl;
    m_vn_vt_unsupported_printed = true;
  }
  std::istringstream zstream(vertex_spec[0]);
  if (!(zstream >> zidx)) {
    std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
              << ", z vertex index must be an integer." << std::endl;
    m_successful_parse = false;
    return;
  }
  if (zidx > (int)tri_mesh.nv()) {
    std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
              << ", z vertex index greater than number of vertices."
              << std::endl;
    m_successful_parse = false;
    return;
  }
  if (zidx < 0)
    if (zidx + ((int)tri_mesh.nv()) < 0) {
      std::cout << "OBJPARSER: Invalid face command on line " << m_line_num
                << ", relative z vertex index out of bounds." << std::endl;
      m_successful_parse = false;
      return;
    }

  //    TriangularMesh::vertex_handle vhx(xidx-1);
  //    TriangularMesh::vertex_handle vhy(yidx-1);
  //    TriangularMesh::vertex_handle vhz(zidx-1);
  //
  //    // Generate three edges for the face
  //    tri_mesh.addEdge(vhx,vhy);
  //    tri_mesh.addEdge(vhy,vhz);
  //    tri_mesh.addEdge(vhz,vhx);
  //
  //    // Generate a triangular face
  //    tri_mesh.addFace(vhx,vhy,vhz);

  // Generate a triangular face
  tri_mesh.addFace(xidx - 1, yidx - 1, zidx - 1);
}

// Code copied from
// http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
void ObjParser::tokenize(const std::string& str,
                         std::vector<std::string>& tokens,
                         const std::string& delimiters) {
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

bool CoffParser::loadParticles(const std::string& coff_file_name,
                               std::vector<Vec3x>& particles) {
  std::ifstream ifs(coff_file_name.c_str());
  if (!ifs.is_open()) {
    std::cout << "COFFPARSER: Failed to open file " << coff_file_name
              << std::endl;
    return false;
  }

  std::string magic;
  ifs >> magic;
  if (magic.compare("COFF") != 0) {
    ifs.close();
    return false;
  }
  int N, b0, b1;
  ifs >> N >> b0 >> b1;
  particles.resize(N);

  for (int i = 0; i < N && !ifs.eof(); ++i) {
    ifs >> particles[i][0] >> particles[i][1] >> particles[i][2];
  }
  ifs.close();
  return true;
}
template <typename S, typename T>
bool ArrayWriter<S, T>::writeArray3(
    const std::string& array_file_name,
    const bridson::Array3<S, bridson::Array1<S> >& arr) {
#ifdef _DEBUG
  std::ofstream ofs(array_file_name.c_str(), std::ios::binary);
  ofs.write(reinterpret_cast<const char*>(&arr.ni), sizeof(arr.ni));
  ofs.write(reinterpret_cast<const char*>(&arr.nj), sizeof(arr.nj));
  ofs.write(reinterpret_cast<const char*>(&arr.nk), sizeof(arr.nk));
  int N = arr.a.size();
  ofs.write(reinterpret_cast<const char*>(&N), sizeof(N));

  std::vector<T> conv_a(N);
  for (int i = 0; i < N; ++i) conv_a[i] = (T)arr.a[i];
  ofs.write(reinterpret_cast<const char*>(&conv_a[0]), N * sizeof(T));
  ofs.close();
#endif
  return true;
}
template <typename S, typename T>
bool ArrayWriter<S, T>::writeArray3(const std::string& array_file_name,
                                    const std::vector<S>& arr, int ni, int nj,
                                    int nk) {
#ifdef _DEBUG
  std::ofstream ofs(array_file_name.c_str(), std::ios::binary);
  ofs.write(reinterpret_cast<const char*>(&ni), sizeof(ni));
  ofs.write(reinterpret_cast<const char*>(&nj), sizeof(nj));
  ofs.write(reinterpret_cast<const char*>(&nk), sizeof(nk));
  int N = arr.size();
  ofs.write(reinterpret_cast<const char*>(&N), sizeof(N));

  std::vector<T> conv_a(N);
  for (int i = 0; i < N; ++i) conv_a[i] = (T)arr[i];
  ofs.write(reinterpret_cast<const char*>(&conv_a[0]), N * sizeof(T));
  ofs.close();
#endif
  return true;
}

template bool ArrayWriter<Scalar>::writeArray3(
    const std::string& array_file_name, const std::vector<Scalar>& arr, int ni,
    int nj, int nk);
template bool ArrayWriter<Scalar>::writeArray3(
    const std::string& array_file_name,
    const bridson::Array3<Scalar, bridson::Array1<Scalar> >& arr);
template bool ArrayWriter<char, Scalar>::writeArray3(
    const std::string& array_file_name, const std::vector<char>& arr, int ni,
    int nj, int nk);
template bool ArrayWriter<char, Scalar>::writeArray3(
    const std::string& array_file_name,
    const bridson::Array3<char, bridson::Array1<char> >& arr);

}  // namespace strandsim
