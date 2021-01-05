/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "XMLReader.hh"

#include "../../../StrandSim/Core/CollisionParameters.hh"
#include "../../../StrandSim/Dynamic/DistanceFieldObject.hh"
#include "../../../StrandSim/Dynamic/ImplicitStepper.hh"
#include "../../../StrandSim/Dynamic/LiquidSimulator.hh"
#include "../../../StrandSim/Forces/FluidDragForce.hh"
#include "../../../StrandSim/Forces/FluidPressureForce.hh"
#include "../../../StrandSim/Utils/StringUtilities.hh"
#include "../../../StrandSim/Utils/ThreadUtils.hh"
#include "../Render/FluidsRenderer.hh"

//#include <boost/filesystem.hpp>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <fstream>

using namespace strandsim;

XMLReader::XMLReader() : ProblemStepper("XMLReader", "read from XML files") {}

XMLReader::~XMLReader() {}

void XMLReader::loadXMLFile(const std::string& filename,
                            std::vector<char>& xmlchars,
                            rapidxml::xml_document<>& doc) {
  // Attempt to read the text from the user-specified xml file
  std::string filecontents;
  if (!loadTextFileIntoString(filename, filecontents)) {
    std::cerr << "ERROR IN TWODSCENEXMLPARSER: XML scene file " << filename
              << ". Failed to read file." << std::endl;
    exit(1);
  }

  m_filename = filename;

  // Copy string into an array of characters for the xml parser
  for (int i = 0; i < (int)filecontents.size(); ++i)
    xmlchars.push_back(filecontents[i]);
  xmlchars.push_back('\0');

  // Initialize the xml parser with the character vector
  doc.parse<0>(&xmlchars[0]);
}

bool XMLReader::loadTextFileIntoString(const std::string& filename,
                                       std::string& filecontents) {
  // Attempt to open the text file for reading
  std::ifstream textfile(filename.c_str(), std::ifstream::in);
  if (!textfile) return false;

  // Read the entire file into a single string
  std::string line;
  while (getline(textfile, line)) filecontents.append(line);

  textfile.close();

  return true;
}

void XMLReader::loadCollisionFree(rapidxml::xml_node<>* node) {
  m_collision_free.clear();
  int forcenum = 0;

  for (rapidxml::xml_node<>* nd = node->first_node("collisionfree"); nd;
       nd = nd->next_sibling("collisionfree")) {
    std::vector<int> cf_arr;

    for (rapidxml::xml_node<>* subnd = nd->first_node("p"); subnd;
         subnd = subnd->next_sibling("p")) {
      int i = -1;
      if (subnd->first_attribute("i")) {
        std::string attribute(subnd->first_attribute("i")->value());
        if (!stringutils::extractFromString(attribute, i)) {
          std::cerr
              << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " Failed to parse value of i attribute for collisionfree "
              << forcenum << ". Value must be integer. Exiting." << std::endl;
          exit(1);
        }
      }

      if (i < 0 || i >= (int)m_global_to_local.size()) continue;
      cf_arr.push_back(i);
    }

    for (int i : cf_arr) {
      const std::pair<int, int>& vidx = m_global_to_local[i];
      for (int j : cf_arr) {
        if (i == j) continue;

        const std::pair<int, int>& nvidx = m_global_to_local[j];
        m_collision_free[vidx].insert(nvidx);
        m_collision_free[nvidx].insert(vidx);
      }
    }

    ++forcenum;
  }
}

int XMLReader::LoadOptions(const char* filename) {
  std::string file_name = std::string(filename);

  size_t start = file_name.find_last_of('/');
  if (start == std::string::npos)
    start = 0;
  else
    start++;

#ifdef WIN32
  size_t wstart = file_name.find_last_of('\\');
  if (wstart == std::string::npos)
    wstart = 0;
  else
    wstart++;
  start = std::max(start, wstart);
#endif

  size_t end = file_name.find_last_of('.');
  if (end == std::string::npos) end = file_name.size();

  size_t count = end - start;

  m_problemName = file_name.substr(start, count);

  loadXMLFile(file_name, m_xmlchars, m_doc);

  m_scene_node = m_doc.first_node("scene");
  if (m_scene_node == NULL) {
    std::cerr << "ERROR IN XMLSCENEPARSER: Failed to parse xml scene file. "
                 "Failed to locate root <scene> node. Exiting."
              << std::endl;
    exit(1);
  }

  loadIntegrator(m_scene_node, m_dt);

  std::cout << "# step_size: \n" << m_dt << std::endl;

  return 0;
}

void XMLReader::setupStrands() {
  int mg_part;
  int mg_df;
  loadParticles(m_scene_node, mg_part);
  loadBucketInfo(m_scene_node);
  loadSolidObject(m_scene_node, mg_df, m_dt);

  loadStrandParameters(m_scene_node, m_dt);

  initGroups(std::max(mg_part, mg_df) + 1);
  loadHairs(m_scene_node);
  loadHairPose(m_scene_node);
  initParticleFaceMapping();

  loadCollisionFree(m_scene_node);

  loadSimpleGravityForces(m_scene_node);
  loadLiquid(m_scene_node, m_dt);

  // init groups

  loadScripts(m_scene_node);
}

void XMLReader::setupAfterInit(int& current_frame, int& current_check_point) {
  // load checkpoint
  loadCheckpoint(m_scene_node, current_frame, current_check_point);
}

void XMLReader::loadSolidObject(rapidxml::xml_node<>* node, int& maxgroup,
                                const double& dt) {
  maxgroup = 0;

  int num_components = 1;

  rapidxml::xml_node<>* nd = node->first_node("liquidinfo");
  if (nd) {
    loadParam(nd, "numComponents", num_components);
  }

  m_num_components = num_components;

  for (rapidxml::xml_node<>* subnd = node->first_node("distancefield"); subnd;
       subnd = subnd->next_sibling("distancefield")) {
    DISTANCE_FIELD_USAGE dfu = DFU_SOLID;
    if (subnd->first_attribute("usage")) {
      std::string handlertype(subnd->first_attribute("usage")->value());
      if (handlertype == "solid")
        dfu = DFU_SOLID;
      else if (handlertype == "source")
        dfu = DFU_SOURCE;
      else if (handlertype == "terminator")
        dfu = DFU_TERMINATOR;
      else {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of type attribute for "
                     "distancefield parameters. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    DISTANCE_FIELD_TYPE bt = DFT_COUNT;
    if (subnd->first_attribute("type")) {
      std::string handlertype(subnd->first_attribute("type")->value());
      if (handlertype == "sphere")
        bt = DFT_SPHERE;
      else if (handlertype == "box")
        bt = DFT_BOX;
      else if (handlertype == "capsule")
        bt = DFT_CAPSULE;
      else if (handlertype == "cylinder")
        bt = DFT_CYLINDER;
      else if (handlertype == "file")
        bt = DFT_FILE;
      else if (handlertype == "sequence")
        bt = DFT_SEQUENCE;
      else {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of type attribute for "
                     "distancefield parameters. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    int color_index = 0;
    if (subnd->first_attribute("colorIndex")) {
      std::string attribute(subnd->first_attribute("colorIndex")->value());
      if (!stringutils::extractFromString(attribute, color_index)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of colorIndex attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    int group = 0;
    if (subnd->first_attribute("group")) {
      std::string attribute(subnd->first_attribute("group")->value());
      if (!stringutils::extractFromString(attribute, group)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of group attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    maxgroup = std::max(maxgroup, group);

    int params_index = 0;
    if (subnd->first_attribute("params")) {
      std::string attribute(subnd->first_attribute("params")->value());
      if (!stringutils::extractFromString(attribute, params_index)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of center(0) attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    bool inverted = false;
    if (subnd->first_attribute("inverted")) {
      std::string attribute(subnd->first_attribute("inverted")->value());
      if (!stringutils::extractFromString(attribute, inverted)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of inverted attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    bool surfflow = false;
    if (subnd->first_attribute("surfflow")) {
      std::string attribute(subnd->first_attribute("surfflow")->value());
      if (!stringutils::extractFromString(attribute, surfflow)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of surfflow attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    bool clamped = true;
    if (subnd->first_attribute("clamped")) {
      std::string attribute(subnd->first_attribute("clamped")->value());
      if (!stringutils::extractFromString(attribute, clamped)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of clamped attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    Scalar flow_height = 0.;
    if (subnd->first_attribute("flowheight")) {
      std::string attribute(subnd->first_attribute("flowheight")->value());
      if (!stringutils::extractFromString(attribute, flow_height)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of flowheight attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    std::vector<DistanceFieldObject::EMIT_INFO> emits;

    if (dfu == DFU_SOURCE) {
      for (rapidxml::xml_node<>* subsubnd = subnd->first_node("emit"); subsubnd;
           subsubnd = subsubnd->next_sibling("emit")) {
        DistanceFieldObject::EMIT_INFO info;
        info.enabled = false;
        info.maxvol = std::numeric_limits<Scalar>::max();

        if (subsubnd->first_attribute("start")) {
          std::string attribute(subsubnd->first_attribute("start")->value());
          if (!stringutils::extractFromString(attribute, info.start)) {
            std::cerr
                << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of start attribute for "
                   "distancefield parameters. Value must be numeric. Exiting."
                << std::endl;
            exit(1);
          }
        }

        if (subsubnd->first_attribute("end")) {
          std::string attribute(subsubnd->first_attribute("end")->value());
          if (!stringutils::extractFromString(attribute, info.end)) {
            std::cerr
                << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of end attribute for distancefield "
                   "parameters. Value must be numeric. Exiting."
                << std::endl;
            exit(1);
          }
        }

        if (subsubnd->first_attribute("evx")) {
          std::string attribute(subsubnd->first_attribute("evx")->value());
          if (!stringutils::extractFromString(attribute, info.emit_vel(0))) {
            std::cerr
                << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of ev(0) attribute for "
                   "distancefield parameters. Value must be numeric. Exiting."
                << std::endl;
            exit(1);
          }
        }

        if (subsubnd->first_attribute("evy")) {
          std::string attribute(subsubnd->first_attribute("evy")->value());
          if (!stringutils::extractFromString(attribute, info.emit_vel(1))) {
            std::cerr
                << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of ev(1) attribute for "
                   "distancefield parameters. Value must be numeric. Exiting."
                << std::endl;
            exit(1);
          }
        }

        if (subsubnd->first_attribute("evz")) {
          std::string attribute(subsubnd->first_attribute("evz")->value());
          if (!stringutils::extractFromString(attribute, info.emit_vel(2))) {
            std::cerr
                << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of ev(2) attribute for "
                   "distancefield parameters. Value must be numeric. Exiting."
                << std::endl;
            exit(1);
          }
        }

        if (subsubnd->first_attribute("maxvol")) {
          std::string attribute(subsubnd->first_attribute("maxvol")->value());
          if (!stringutils::extractFromString(attribute, info.maxvol)) {
            std::cerr
                << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of maxvol attribute for "
                   "distancefield parameters. Value must be numeric. Exiting."
                << std::endl;
            exit(1);
          }
        }

        emits.push_back(info);
      }

      if (emits.size() == 0) {
        Vec3x eject_vel = Vec3x::Zero();
        if (subnd->first_attribute("vx")) {
          std::string attribute(subnd->first_attribute("vx")->value());
          if (!stringutils::extractFromString(attribute, eject_vel(0))) {
            std::cerr << outputmod::startred
                      << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                      << " Failed to parse end attribute for vx parameters. "
                         "Value must be numeric. Exiting."
                      << std::endl;
            exit(1);
          }
        }

        if (subnd->first_attribute("vy")) {
          std::string attribute(subnd->first_attribute("vy")->value());
          if (!stringutils::extractFromString(attribute, eject_vel(1))) {
            std::cerr << outputmod::startred
                      << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                      << " Failed to parse end attribute for vy parameters. "
                         "Value must be numeric. Exiting."
                      << std::endl;
            exit(1);
          }
        }

        if (subnd->first_attribute("vz")) {
          std::string attribute(subnd->first_attribute("vz")->value());
          if (!stringutils::extractFromString(attribute, eject_vel(2))) {
            std::cerr << outputmod::startred
                      << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                      << " Failed to parse end attribute for vz parameters. "
                         "Value must be numeric. Exiting."
                      << std::endl;
            exit(1);
          }
        }
        DistanceFieldObject::EMIT_INFO emit = {
            eject_vel, 0.0, 0.0, std::numeric_limits<Scalar>::max(), true};
        emits.push_back(emit);
      }
    }

    if (bt == DFT_BOX || bt == DFT_SPHERE || bt == DFT_CAPSULE ||
        bt == DFT_CYLINDER || bt == DFT_FILE || bt == DFT_SEQUENCE) {
      Vec3x center;
      Vec3x raxis = Vec3x(0, 1, 0);
      double rangle = 0.0;

      if (subnd->first_attribute("cx")) {
        std::string attribute(subnd->first_attribute("cx")->value());
        if (!stringutils::extractFromString(attribute, center(0))) {
          std::cerr
              << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " Failed to parse value of center(0) attribute for "
                 "distancefield parameters. Value must be numeric. Exiting."
              << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("cy")) {
        std::string attribute(subnd->first_attribute("cy")->value());
        if (!stringutils::extractFromString(attribute, center(1))) {
          std::cerr
              << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " Failed to parse value of center(1) attribute for "
                 "distancefield parameters. Value must be numeric. Exiting."
              << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("cz")) {
        std::string attribute(subnd->first_attribute("cz")->value());
        if (!stringutils::extractFromString(attribute, center(2))) {
          std::cerr
              << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " Failed to parse value of center(2) attribute for "
                 "distancefield parameters. Value must be numeric. Exiting."
              << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("rx")) {
        std::string attribute(subnd->first_attribute("rx")->value());
        if (!stringutils::extractFromString(attribute, raxis(0))) {
          std::cerr
              << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " Failed to parse value of center(0) attribute for "
                 "distancefield parameters. Value must be numeric. Exiting."
              << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("ry")) {
        std::string attribute(subnd->first_attribute("ry")->value());
        if (!stringutils::extractFromString(attribute, raxis(1))) {
          std::cerr
              << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " Failed to parse value of center(1) attribute for "
                 "distancefield parameters. Value must be numeric. Exiting."
              << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("rz")) {
        std::string attribute(subnd->first_attribute("rz")->value());
        if (!stringutils::extractFromString(attribute, raxis(2))) {
          std::cerr
              << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " Failed to parse value of center(2) attribute for "
                 "distancefield parameters. Value must be numeric. Exiting."
              << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("rw")) {
        std::string attribute(subnd->first_attribute("rw")->value());
        if (!stringutils::extractFromString(attribute, rangle)) {
          std::cerr
              << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " Failed to parse value of center(2) attribute for "
                 "distancefield parameters. Value must be numeric. Exiting."
              << std::endl;
          exit(1);
        }
      }

      Vec4x parameter(0, 0, 0, 0);
      Scalar dx_seq = m_dx;

      switch (bt) {
        case DFT_SPHERE:
          if (subnd->first_attribute("radius")) {
            std::string attribute(subnd->first_attribute("radius")->value());
            if (!stringutils::extractFromString(attribute, parameter(0))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of radius attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          break;
        case DFT_BOX:
          if (subnd->first_attribute("ex")) {
            std::string attribute(subnd->first_attribute("ex")->value());
            if (!stringutils::extractFromString(attribute, parameter(0))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of ex attribute for distancefield "
                     "parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }

          if (subnd->first_attribute("ey")) {
            std::string attribute(subnd->first_attribute("ey")->value());
            if (!stringutils::extractFromString(attribute, parameter(1))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of ey attribute for distancefield "
                     "parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }

          if (subnd->first_attribute("ez")) {
            std::string attribute(subnd->first_attribute("ez")->value());
            if (!stringutils::extractFromString(attribute, parameter(2))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of ez attribute for distancefield "
                     "parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }

          if (subnd->first_attribute("radius")) {
            std::string attribute(subnd->first_attribute("radius")->value());
            if (!stringutils::extractFromString(attribute, parameter(3))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of radius attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          break;
        case DFT_CAPSULE:
          if (subnd->first_attribute("radius")) {
            std::string attribute(subnd->first_attribute("radius")->value());
            if (!stringutils::extractFromString(attribute, parameter(0))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of radius attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          if (subnd->first_attribute("halflength")) {
            std::string attribute(
                subnd->first_attribute("halflength")->value());
            if (!stringutils::extractFromString(attribute, parameter(1))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of halflength attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          break;
        case DFT_CYLINDER:
          if (subnd->first_attribute("radius")) {
            std::string attribute(subnd->first_attribute("radius")->value());
            if (!stringutils::extractFromString(attribute, parameter(0))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of radius attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          if (subnd->first_attribute("corner")) {
            std::string attribute(subnd->first_attribute("corner")->value());
            if (!stringutils::extractFromString(attribute, parameter(1))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of corner attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          if (subnd->first_attribute("halflength")) {
            std::string attribute(
                subnd->first_attribute("halflength")->value());
            if (!stringutils::extractFromString(attribute, parameter(2))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of halflength attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          break;
        case DFT_FILE:
          parameter(0) = 1.0;
          if (subnd->first_attribute("scale")) {
            std::string attribute(subnd->first_attribute("scale")->value());
            if (!stringutils::extractFromString(attribute, parameter(0))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of scale attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          if (subnd->first_attribute("dx")) {
            std::string attribute(subnd->first_attribute("dx")->value());
            if (!stringutils::extractFromString(attribute, dx_seq)) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of dx attribute for distancefield "
                     "parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          break;
        case DFT_SEQUENCE:
          parameter = Vec4x(0, 1, 0, 1);
          if (subnd->first_attribute("startTime")) {
            std::string attribute(subnd->first_attribute("startTime")->value());
            if (!stringutils::extractFromString(attribute, parameter(0))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of startTime attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          if (subnd->first_attribute("endTime")) {
            std::string attribute(subnd->first_attribute("endTime")->value());
            if (!stringutils::extractFromString(attribute, parameter(1))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of endTime attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          if (subnd->first_attribute("startFrame")) {
            std::string attribute(
                subnd->first_attribute("startFrame")->value());
            if (!stringutils::extractFromString(attribute, parameter(2))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of startFrame attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          if (subnd->first_attribute("endFrame")) {
            std::string attribute(subnd->first_attribute("endFrame")->value());
            if (!stringutils::extractFromString(attribute, parameter(3))) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of endFrame attribute for "
                     "distancefield parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          if (subnd->first_attribute("dx")) {
            std::string attribute(subnd->first_attribute("dx")->value());
            if (!stringutils::extractFromString(attribute, dx_seq)) {
              std::cerr
                  << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of dx attribute for distancefield "
                     "parameters. Value must be numeric. Exiting."
                  << std::endl;
              exit(1);
            }
          }
          break;
        default:
          break;
      }

      if (bt == DFT_SEQUENCE) {
        std::string filename;

        if (subnd->first_attribute("base")) {
          filename = std::string(subnd->first_attribute("base")->value());
        }

        if (!filename.empty()) {
          if (dfu == DFU_SOLID) {
            m_fields.push_back(DistanceFieldObject(
                center, parameter, bt, dfu, raxis, rangle, group, color_index,
                dt, dx_seq, emits, inverted, false, surfflow, filename, ""));
            m_solid_flow_height.push_back(flow_height);
          } else if (dfu == DFU_SOURCE) {
            m_sources.push_back(DistanceFieldObject(
                center, parameter, bt, dfu, raxis, rangle, group, color_index,
                dt, dx_seq, emits, inverted, clamped, surfflow, filename, ""));
          } else if (dfu == DFU_TERMINATOR) {
            m_terminators.push_back(DistanceFieldObject(
                center, parameter, bt, dfu, raxis, rangle, group, color_index,
                dt, dx_seq, emits, inverted, false, surfflow, filename, ""));
          }
        } else {
          continue;
        }

      } else if (bt == DFT_FILE) {
        std::string filename;

        if (subnd->first_attribute("filename")) {
          filename = std::string(subnd->first_attribute("filename")->value());
        }

        bool cached = false;
        if (subnd->first_attribute("cached")) {
          std::string attribute(subnd->first_attribute("cached")->value());
          if (!stringutils::extractFromString(attribute, cached)) {
            std::cerr
                << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of cached attribute for "
                   "distancefield parameters. Value must be numeric. Exiting."
                << std::endl;
            exit(1);
          }
        }

        std::string cachename;

        if (cached) {
          cachename = filename + ".cache";
          if (subnd->first_attribute("cachename")) {
            cachename =
                std::string(subnd->first_attribute("cachename")->value());
          }
        }

        if (!filename.empty()) {
          if (dfu == DFU_SOLID) {
            m_fields.push_back(DistanceFieldObject(
                center, parameter, bt, dfu, raxis, rangle, group, color_index,
                dt, dx_seq, emits, inverted, false, surfflow, filename,
                cachename));
            m_solid_flow_height.push_back(flow_height);
          } else if (dfu == DFU_SOURCE) {
            m_sources.push_back(DistanceFieldObject(
                center, parameter, bt, dfu, raxis, rangle, group, color_index,
                dt, dx_seq, emits, inverted, clamped, surfflow, filename,
                cachename));
          } else if (dfu == DFU_TERMINATOR) {
            m_terminators.push_back(DistanceFieldObject(
                center, parameter, bt, dfu, raxis, rangle, group, color_index,
                dt, dx_seq, emits, inverted, false, surfflow, filename,
                cachename));
          }
        } else {
          continue;
        }

      } else {
        if (dfu == DFU_SOLID) {
          m_fields.push_back(DistanceFieldObject(
              center, parameter, bt, dfu, raxis, rangle, group, color_index, dt,
              dx_seq, emits, inverted, false, surfflow));
          m_solid_flow_height.push_back(flow_height);
        } else if (dfu == DFU_SOURCE) {
          m_sources.push_back(DistanceFieldObject(
              center, parameter, bt, dfu, raxis, rangle, group, color_index, dt,
              dx_seq, emits, inverted, clamped, surfflow));
        } else {
          m_terminators.push_back(DistanceFieldObject(
              center, parameter, bt, dfu, raxis, rangle, group, color_index, dt,
              dx_seq, emits, inverted, false, surfflow));
        }
      }

      if (dfu == DFU_SOLID) {
        DistanceFieldObject& dfo = m_fields[m_fields.size() - 1];

        m_meshScripting_controllers.push_back(dfo.mesh_controller);

        dfo.mesh_controller->setDefaultFrictionCoefficient(
            m_simulation_params.m_hairMeshFrictionCoefficient);

        auto currentMesh = dfo.mesh_controller->getCurrentMesh();

        transformTriangleObject(*currentMesh, dfo.rot, Vec3x::Zero(),
                                dfo.center);

        dfo.resetDisplacement();

        strandsim::TriangleMeshRenderer* mesh_renderer =
            new strandsim::TriangleMeshRenderer(
                *(dfo.mesh_controller->getCurrentMesh()), num_components);
        m_mesh_renderers.push_back(mesh_renderer);
      } else if (dfu == DFU_SOURCE) {
        DistanceFieldObject& dfo = m_sources[m_sources.size() - 1];

        auto currentMesh = dfo.mesh_controller->getCurrentMesh();

        transformTriangleObject(*currentMesh, dfo.rot, Vec3x::Zero(),
                                dfo.center);

        dfo.resetDisplacement();
      } else if (dfu == DFU_TERMINATOR) {
        DistanceFieldObject& dfo = m_terminators[m_terminators.size() - 1];

        auto currentMesh = dfo.mesh_controller->getCurrentMesh();

        transformTriangleObject(*currentMesh, dfo.rot, Vec3x::Zero(),
                                dfo.center);

        dfo.resetDisplacement();

        strandsim::TriangleMeshRenderer* mesh_renderer =
            new strandsim::TriangleMeshRenderer(
                *(dfo.mesh_controller->getCurrentMesh()), num_components);
        m_mesh_renderers.push_back(mesh_renderer);
      }
    }
  }
}

void XMLReader::loadSimpleGravityForces(rapidxml::xml_node<>* node) {
  assert(node != NULL);

  // Load each constant force
  rapidxml::xml_node<>* nd = node->first_node("simplegravity");

  m_gravity.setZero();

  if (nd) {
    // Extract the x component of the force
    if (nd->first_attribute("fx")) {
      std::string attribute(nd->first_attribute("fx")->value());
      if (!stringutils::extractFromString(attribute, m_gravity.x())) {
        exit(1);
      }
    }

    // Extract the y component of the force
    if (nd->first_attribute("fy")) {
      std::string attribute(nd->first_attribute("fy")->value());
      if (!stringutils::extractFromString(attribute, m_gravity.y())) {
        exit(1);
      }
    }

    // Extract the z component of the force
    if (nd->first_attribute("fz")) {
      std::string attribute(nd->first_attribute("fz")->value());
      if (!stringutils::extractFromString(attribute, m_gravity.z())) {
        exit(1);
      }
    }
  }

  GravitationForce::setGravity(m_gravity);
}

void XMLReader::loadBucketInfo(rapidxml::xml_node<>* node) {
  rapidxml::xml_node<>* nd = node->first_node("bucketinfo");
  m_bucket_size = 5.0;
  m_num_cells = 5;
  m_dx = 1.0;

  if (!nd) return;

  if (nd) {
    if (nd->first_attribute("size")) {
      std::string attribute(nd->first_attribute("size")->value());
      if (!stringutils::extractFromString(attribute, m_bucket_size)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of size for bucketinfo. Value "
                     "must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("numcells")) {
      std::string attribute(nd->first_attribute("numcells")->value());
      if (!stringutils::extractFromString(attribute, m_num_cells)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of numcells for bucketinfo. Value "
                     "must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }
  }

  m_dx = m_bucket_size / (Scalar)m_num_cells;
}

void XMLReader::PrintAdditionalSettings(const std::string& szfn_dir) {}

template <typename T>
void XMLReader::loadParam(rapidxml::xml_node<>* nd, const char* name,
                          Eigen::Matrix<T, Eigen::Dynamic, 1>& param,
                          int num_components) {
  rapidxml::xml_node<>* subnd = NULL;
  if ((subnd = nd->first_node(name))) {
    if (subnd->first_attribute("value")) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, param(0))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of " << name << " attribute for "
                  << nd->name() << ". Exiting." << std::endl;
        exit(1);
      }
      param = Eigen::Matrix<T, Eigen::Dynamic, 1>::Constant(num_components,
                                                            param(0));
    } else {
      int i = 0;
      param = Eigen::Matrix<T, Eigen::Dynamic, 1>::Constant(num_components,
                                                            param(0));

      for (rapidxml::xml_node<>* subsubnd = subnd->first_node("component");
           subsubnd; subsubnd = subsubnd->next_sibling("component")) {
        if (subsubnd->first_attribute("value")) {
          std::string attribute(subsubnd->first_attribute("value")->value());
          if (!stringutils::extractFromString(attribute, param(i))) {
            std::cerr << outputmod::startred
                      << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                      << " Failed to parse " << i << "-th value of " << name
                      << " attribute for " << nd->name() << ". Exiting."
                      << std::endl;
            exit(1);
          }
        }
        ++i;
      }
    }
  } else {
    param =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Constant(num_components, param(0));
  }
}

void XMLReader::loadLiquid(rapidxml::xml_node<>* node, const double& dt) {
  LiquidInfo info;

  info.num_components = 1;

  info.liquid_density = VecXx::Constant(1, 1.0);  // g/cm^3
  info.surf_tension_coeff = VecXx::Constant(1, 72.8);
  info.liquid_bulk_modulus = VecXx::Constant(1, 1.09e6);
  info.liquid_shear_modulus = VecXx::Constant(1, 2.9e3);
  info.plastic_weaken_strain = VecXx::Constant(1, 2.175e2);
  info.plastic_yield_stress = VecXx::Constant(1, 3.19e2);
  info.flow_consistency_index = VecXx::Constant(1, 2.72e2);
  info.flow_behavior_index = VecXx::Constant(1, 0.22);
  info.plastic_relax_coeff = VecXx::Constant(1, 0.35);
  info.helmholtz_first_minima = VecXx::Constant(1, 0.5);
  info.helmholtz_second_minima = VecXx::Constant(1, 0.5);

  info.solid_shell_thickness = 1.0;
  info.correction_step = 1;
  info.correction_multiplier = 2.0;
  info.correction_strength = 10.0;
  info.affine_stretch_damping = 0.0;
  info.affine_rotate_damping = 0.0;
  info.velocity_damping = 0.0;
  info.particle_cell_multiplier = 0.5;
  info.liquid_boundary_friction = 1.0;
  info.surf_tension_smoothing_step = 7;
  info.iteration_print_step = 0;
  info.elasto_capture_rate = 1.0;
  info.theta_criterion = 0.01;
  info.shear_pcg_criterion = 1e-4;
  info.pressure_pcg_criterion = 1e-8;
  info.pcg_max_iters = 1000;
  info.typical_flow_thickness = 0.05;
  info.levelset_young_modulus = 0.0;
  info.mesh_flow_epsilon = 0.002;
  info.mesh_flow_slip_length = 0.0002;
  info.mesh_flow_max_height = 0.5;
  info.mesh_flow_critical_height = 0.25;
  info.liquid_shear_damping = 0.001;
  info.rest_contact_angle = 60.0 / 180.0 * M_PI;
  info.chemical_diffusivity = 1.0;
  info.geometric_diffusivity = 1.0;
  info.signed_distance_multiplier = 2;
  info.geometric_drag_insulation = 0.0;

  info.use_surf_tension = false;
  info.use_implicit_elasticity = true;
  info.use_liquid_capture = true;
  info.use_liquid_drag = true;
  info.use_varying_volume_fraction = true;
  info.solid_projection = true;
  info.solid_adhesion = true;
  info.solve_viscosity = true;
  info.solve_color_diffusion = false;
  info.use_implicit_pressure = true;
  info.use_constant_drag = false;

  rapidxml::xml_node<>* nd = node->first_node("liquidinfo");
  if (nd) {
    rapidxml::xml_node<>* subnd = NULL;

    loadParam(nd, "numComponents", info.num_components);
    loadParam(nd, "liquidDensity", info.liquid_density, info.num_components);
    loadParam(nd, "surfTensionCoeff", info.surf_tension_coeff,
              info.num_components);
    loadParam(nd, "plasticWeakenStrain", info.plastic_weaken_strain,
              info.num_components);
    loadParam(nd, "liquidBulkModulus", info.liquid_bulk_modulus,
              info.num_components);
    loadParam(nd, "flowBehaviorIndex", info.flow_behavior_index,
              info.num_components);
    loadParam(nd, "flowConsistencyIndex", info.flow_consistency_index,
              info.num_components);
    loadParam(nd, "plasticYieldStress", info.plastic_yield_stress,
              info.num_components);
    loadParam(nd, "liquidShearModulus", info.liquid_shear_modulus,
              info.num_components);
    loadParam(nd, "plasticRelaxCoeff", info.plastic_relax_coeff,
              info.num_components);
    loadParam(nd, "helmholtzFirstMinima", info.helmholtz_first_minima,
              info.num_components);
    loadParam(nd, "helmholtzSecondMinima", info.helmholtz_second_minima,
              info.num_components);

    loadParam(nd, "pcgMaxIters", info.pcg_max_iters);
    loadParam(nd, "chemicalDiffusivity", info.chemical_diffusivity);
    loadParam(nd, "geometricDiffusivity", info.geometric_diffusivity);
    loadParam(
        nd, "restContactAngle", info.rest_contact_angle,
        std::function<double(const double&)>(
            [](const double& val) -> double { return val * M_PI / 180.0; }));
    loadParam(nd, "liquidShearDamping", info.liquid_shear_damping);
    loadParam(nd, "signedDistanceMultiplier", info.signed_distance_multiplier);
    loadParam(nd, "geometricDragInsulation", info.geometric_drag_insulation);
    loadParam(nd, "shearPcgCriterion", info.shear_pcg_criterion);
    loadParam(nd, "pressurePcgCriterion", info.pressure_pcg_criterion);
    loadParam(nd, "solveColorDiffusion", info.solve_color_diffusion);
    loadParam(nd, "solveViscosity", info.solve_viscosity);
    loadParam(nd, "liquidCapture", info.use_liquid_capture);
    loadParam(nd, "liquidDrag", info.use_liquid_drag);
    loadParam(nd, "varyingVolumeFraction", info.use_varying_volume_fraction);
    loadParam(nd, "constantDrag", info.use_constant_drag);
    loadParam(nd, "solidShellThickness", info.solid_shell_thickness);
    loadParam(nd, "meshFlowCriticalHeight", info.mesh_flow_critical_height);
    loadParam(nd, "meshFlowMaxHeight", info.mesh_flow_max_height);
    loadParam(nd, "meshFlowEpsilon", info.mesh_flow_epsilon);
    loadParam(nd, "levelsetYoungModulus", info.levelset_young_modulus);
    loadParam(nd, "typicalFlowThickness", info.typical_flow_thickness);
    loadParam(nd, "liquidBoundaryFriction", info.liquid_boundary_friction);
    loadParam(nd, "thetaCriterion", info.theta_criterion);
    loadParam(nd, "surfTensionSmoothingStep", info.surf_tension_smoothing_step);
    loadParam(nd, "iterationPrintStep", info.iteration_print_step);
    loadParam(nd, "particleCellMultiplier", info.particle_cell_multiplier);
    loadParam(nd, "elastoCaptureRate", info.elasto_capture_rate);
    loadParam(nd, "affineStretchDamping", info.affine_stretch_damping);
    loadParam(nd, "affineRotateDamping", info.affine_rotate_damping);
    loadParam(nd, "velocityDamping", info.velocity_damping);
    loadParam(nd, "solidAdhesion", info.solid_adhesion);
    loadParam(nd, "solidProjection", info.solid_projection);
    loadParam(nd, "surfTension", info.use_surf_tension);
    loadParam(nd, "implicitPressure", info.use_implicit_pressure);
    loadParam(nd, "implicitElasticity", info.use_implicit_elasticity);
    loadParam(nd, "correctionStep", info.correction_step);
    loadParam(nd, "correctionStrength", info.correction_strength);
    loadParam(nd, "correctionMultiplier", info.correction_multiplier);
  }

  std::shared_ptr<FluidScriptingController> controller =
      std::make_shared<LiquidSimulator>(m_fields, m_sources, m_terminators,
                                        m_strands, m_bucket_size, m_num_cells,
                                        dt, info);
  m_fluidScripting_controllers.push_back(controller);

  FluidsRenderer* render = new FluidsRenderer(
      controller, controller->getParticleMutex(), controller->getGridMutex(),
      FluidsRenderer::DrawMode::DBG);
  m_fluids_renderers.push_back(render);

  if (m_fluidScripting_controllers.size() > 0) {
    FluidDragForce::setScriptingController(m_fluidScripting_controllers[0]);
    FluidPressureForce::setScriptingController(m_fluidScripting_controllers[0]);
  }

  for (CollisionParameters& param : m_collision_parameters) {
    param.m_cohesionSigma = info.surf_tension_coeff;
    param.m_cohesionTheta = info.rest_contact_angle;
    param.m_cohesionMaxDist = m_bucket_size / (Scalar)m_num_cells * 0.5;

    param.initializeCohesionTable();
  }

  for (int i = 0; i < (int)m_fields.size(); ++i) {
    m_fields[i].init_mesh_flow(controller, m_solid_flow_height[i]);
  }
}

void XMLReader::initGroups(int num_group) {
  m_groups.resize(num_group);
  m_group_fields.resize(num_group);
  m_group_sources.resize(num_group);
  m_group_pos.resize(num_group, Vec3x::Zero());
  m_group_scale.resize(num_group, Vec3x::Ones());
  m_group_rot.resize(num_group, Eigen::Quaternion<double>::Identity());

  const int num_fields = (int)m_fields.size();
  for (int i = 0; i < num_fields; ++i) {
    m_group_fields[m_fields[i].group].push_back(i);
    if (m_group_fields[m_fields[i].group].size() > 1) {
      std::cout << "Error: More than 1 fields in the same group!" << std::endl;
      exit(-1);
    }
  }

  const int num_sources = (int)m_sources.size();
  for (int i = 0; i < num_sources; ++i) {
    m_group_sources[m_sources[i].group].push_back(i);
    if (m_group_sources[m_sources[i].group].size() > 1) {
      std::cout << "Error: More than 1 source in the same group!" << std::endl;
      exit(-1);
    }

    if (m_group_fields[m_sources[i].group].size() > 0) {
      std::cout << "Error: Source group is overlapping with field group!"
                << std::endl;
      exit(-1);
    }
  }

  // we assume one group has at most one solid object
  for (int i = 0; i < num_group; ++i) {
    if (m_group_fields[i].size() > 0)
      m_group_pos[i] = m_fields[m_group_fields[i][0]].center;
    else if (m_group_sources[i].size() > 0)
      m_group_pos[i] = m_sources[m_group_sources[i][0]].center;
  }

  m_group_prev_pos = m_group_pos;
  m_group_prev_rot = m_group_rot;
  m_group_prev_scale = m_group_scale;
}

void XMLReader::initParticleFaceMapping() {
  const int num_strand = m_local_to_global.size();

  strandsim::for_each(0, num_strand, [&](int sidx) {
    const std::vector<int>& indices = m_local_to_global[sidx];

    int prev_group = -1;
    int prev_closest_face = -1;
    const int num_locals = indices.size();

    for (int i = 0; i < num_locals; ++i) {
      int global_idx = indices[i];

      if (!m_fixed[global_idx]) {
        m_particle_closed_face_index[global_idx] = -1;
        prev_closest_face = -1;
        continue;
      }

      int group = m_particle_groups[global_idx];

      if (prev_closest_face >= 0 && group == prev_group) {
        m_particle_closed_face_index[global_idx] = prev_closest_face;
        continue;
      }

      const std::vector<int>& fields = m_group_fields[group];
      if (fields.size() == 0) {
        m_particle_closed_face_index[global_idx] = -1;
        prev_closest_face = -1;
        continue;
      }

      const DistanceFieldObject& dfo = m_fields[fields[0]];
      Vec3x pos = m_particles[global_idx].segment<3>(0);

      int iFace = dfo.get_closest_face(pos);
      m_particle_closed_face_index[global_idx] = iFace;
      prev_closest_face = iFace;
      prev_group = group;
    }
  });
}

void XMLReader::setupMeshes() {}

void XMLReader::executeHairPose() {
  const double current_time = getTime();
  // first do hair load
  const int num_hairposes = (int)m_hairposes.size();
  for (int i = 0; i < num_hairposes; ++i) {
    if (current_time >= m_hairposes[i].time) {
      if (m_hairposes[i].loaded) continue;

      m_hairposes[i].loaded = true;

      const int num_strands = m_strands.size();

      for (int j = 0; j < num_strands; ++j) {
        m_strands[j]->setCurrentDegreesOfFreedom(m_hairposes[i].dofs[j]);
        m_strands[j]->setFutureDegreesOfFreedom(
            m_strands[j]->getCurrentDegreesOfFreedom());
        m_strands[j]->setCurrentAreaDegreesOfFreedom(m_hairposes[i].areas[j]);
        m_strands[j]->setFutureAreaDegreesOfFreedom(
            m_strands[j]->getCurrentAreaDegreesOfFreedom());

        const int num_verts = m_strands[j]->getNumVertices();
        for (int k = 0; k < num_verts; ++k) {
          if (m_strands[j]->isVertexFreezed(k)) {
            m_strands[j]->dynamics().getScriptingController()->freezeVertices(
                k);
          }

          if (k != num_verts - 1 && m_strands[j]->isThetaFreezed(k)) {
            m_strands[j]->dynamics().getScriptingController()->freezeTheta(k);
          }

          if (m_strands[j]->isVertexGoaled(k)) {
            //                controller->freezeVertices(i);
            Vec3x pos = m_hairposes[i].dofs[j].segment<3>(k * 4);
            m_strands[j]->dynamics().getScriptingController()->setVertexGoal(
                k, pos);
          }

          if (k != num_verts - 1 && m_strands[j]->isThetaGoaled(k)) {
            //                controller->freezeTheta(i);
            m_strands[j]->dynamics().getScriptingController()->setThetaGoal(
                k, m_hairposes[i].dofs[j](k * 4 + 3));
          }
        }
      }
    }
  }
}

bool XMLReader::executeScript(int total_substep_id, const Scalar substep_dt) {
  const int num_scripts = (int)m_scripts.size();
  const double current_time = getTime();

  if (current_time > m_end_time) {
    std::cout << "Simulation Finished!" << std::endl;
    return false;
  }

  executeHairPose();

  for (int i = 0; i < num_scripts; ++i) {
    m_scripts[i].stepScript(substep_dt, current_time);
  }

  for (DistanceFieldObject& obj : m_sources) {
    for (DistanceFieldObject::EMIT_INFO& info : obj.emits) {
      if (current_time >= info.start && current_time < info.end) {
        info.enabled = true;
      } else {
        info.enabled = false;
      }
    }
  }
  // apply the movement
  const int num_groups = (int)m_groups.size();
  for (int i = 0; i < num_groups; ++i) {
    const Eigen::Quaternion<double>& q = m_group_rot[i];
    const Vec3x& t = m_group_pos[i];
    const Vec3x& s = m_group_scale[i];

    const Eigen::Quaternion<double>& q_prev = m_group_prev_rot[i];
    const Vec3x& t_prev = m_group_prev_pos[i];
    const Vec3x& s_prev = m_group_prev_scale[i];

    const Eigen::Quaternion<double> q_diff = q * q_prev.inverse();
    const Vec3x t_diff = t - t_prev;
    const Vec3x s_diff = Vec3x(s.array() / s_prev.array());

    const std::vector<int>& field_indices = m_group_fields[i];
    for (int fidx : field_indices) {
      m_fields[fidx].advance(substep_dt, total_substep_id);
      m_fields[fidx].step_flow_dynamics(substep_dt);
    }

    const std::vector<int>& source_indices = m_group_sources[i];
    for (int fidx : source_indices) {
      m_sources[fidx].advance(substep_dt, total_substep_id);
    }

    const std::vector<std::pair<int, int> >& group = m_groups[i];
    for (const std::pair<int, int>& p : group) {
      const int global_idx = m_local_to_global[p.first][p.second];
      if (m_fixed[global_idx]) {
        if (m_particle_closed_face_index[global_idx] >= 0 &&
            m_group_fields[i].size() > 0) {
          int fidx = m_group_fields[i][0];
          int face = m_particle_closed_face_index[global_idx];

          Vec3x center, translate;
          Eigen::Quaternion<Scalar> rot;
          m_fields[fidx].mesh_controller->getCurrentMesh()->getRigidTransform(
              face, center, translate, rot);
          transformRodRootVtx(*(m_rodDatum[p.first]), rot, center, translate, s,
                              p.second);
        } else {
          transformRodRootVtx(*(m_rodDatum[p.first]), q_diff, t_prev, t_diff,
                              s_diff, p.second);
        }
      }
    }

    m_group_prev_pos[i] = m_group_pos[i];
    m_group_prev_rot[i] = m_group_rot[i];
    m_group_prev_scale[i] = m_group_scale[i];
  }

  return true;
}

void XMLReader::loadParticles(rapidxml::xml_node<>* node, int& maxgroup) {
  int numparticles = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("particle"); nd;
       nd = nd->next_sibling("particle"))
    ++numparticles;

  m_particles.resize(numparticles);
  m_particle_groups.resize(numparticles);
  m_particle_closed_face_index.resize(numparticles);
  m_fixed.resize(numparticles);
  m_surf_height.resize(numparticles);

  maxgroup = 0;

  int particle = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("particle"); nd;
       nd = nd->next_sibling("particle")) {
    // Extract the particle's initial position
    Vec3x pos = Vec3x::Zero();
    if (nd->first_attribute("x")) {
      std::string position(nd->first_attribute("x")->value());
      if (!stringutils::readList(position, ' ', pos)) {
        std::cerr << "Failed to load x, y, and z positions for particle "
                  << particle << std::endl;
        exit(1);
      }
    } else {
      std::cerr
          << "Failed to find x, y, and z position attributes for particle "
          << particle << std::endl;
      exit(1);
    }

    // parse theta
    double theta = 0.0;
    if (nd->first_attribute("theta")) {
      std::string attribute(nd->first_attribute("theta")->value());
      if (!stringutils::extractFromString(attribute, theta)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of fixed attribute for particle "
                  << particle << ". Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    m_particles[particle] = Vec4x(pos(0), pos(1), pos(2), theta);

    // Determine if the particle is fixed
    int fixed = 0;
    if (nd->first_attribute("fixed")) {
      std::string attribute(nd->first_attribute("fixed")->value());
      if (!stringutils::extractFromString(attribute, fixed)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of fixed attribute for particle "
                  << particle << ". Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    m_fixed[particle] = fixed;

    double h = 0.0;
    if (nd->first_attribute("h")) {
      std::string attribute(nd->first_attribute("h")->value());
      if (!stringutils::extractFromString(attribute, h)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse h attribute for particle " << particle
                  << ". Value must be double. Exiting." << std::endl;
        exit(1);
      }
    }
    m_surf_height[particle] = h;

    int group = 0;
    if (nd->first_attribute("group")) {
      std::string attribute(nd->first_attribute("group")->value());
      if (!stringutils::extractFromString(attribute, group)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse group attribute for particle "
                  << particle << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    m_particle_groups[particle] = group;

    m_particle_closed_face_index[particle] = -1;

    maxgroup = std::max(maxgroup, group);

    ++particle;
  }

  // load from obj hairs (used to interacts with Houdini etc.)
  for (rapidxml::xml_node<>* nd = node->first_node("hairobj"); nd;
       nd = nd->next_sibling("hairobj")) {
    std::string szfn;
    if (nd->first_attribute("filename")) {
      szfn = nd->first_attribute("filename")->value();
    }

    if (szfn.empty()) continue;

    std::ifstream ifs(szfn);

    if (ifs.fail()) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to read file: " << szfn << ". Exiting."
                << std::endl;
      exit(1);
    }

    int base_idx = m_particles.size();
    int group = 0;
    if (nd->first_attribute("group")) {
      std::string attribute(nd->first_attribute("group")->value());
      if (!stringutils::extractFromString(attribute, group)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse group attribute for hairobj. Value must "
                     "be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    Vec3x rot_axis = Vec3x::UnitY();
    Scalar rot_angle = 0.0;
    Vec3x translate = Vec3x::Zero();

    if (nd->first_attribute("x")) {
      std::string attribute(nd->first_attribute("x")->value());
      if (!stringutils::extractFromString(attribute, translate(0))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse x attribute for hairobj. Value must be "
                     "scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("y")) {
      std::string attribute(nd->first_attribute("y")->value());
      if (!stringutils::extractFromString(attribute, translate(1))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse y attribute for hairobj. Value must be "
                     "scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("z")) {
      std::string attribute(nd->first_attribute("z")->value());
      if (!stringutils::extractFromString(attribute, translate(2))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse z attribute for hairobj. Value must be "
                     "scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("rx")) {
      std::string attribute(nd->first_attribute("rx")->value());
      if (!stringutils::extractFromString(attribute, rot_axis(0))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse rx attribute for hairobj. Value must be "
                     "scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("ry")) {
      std::string attribute(nd->first_attribute("ry")->value());
      if (!stringutils::extractFromString(attribute, rot_axis(1))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse ry attribute for hairobj. Value must be "
                     "scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("rz")) {
      std::string attribute(nd->first_attribute("rz")->value());
      if (!stringutils::extractFromString(attribute, rot_axis(2))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse rz attribute for hairobj. Value must be "
                     "scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("rw")) {
      std::string attribute(nd->first_attribute("rw")->value());
      if (!stringutils::extractFromString(attribute, rot_angle)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse rw attribute for hairobj. Value must be "
                     "scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    Eigen::Quaternion<Scalar> qrot = Eigen::Quaternion<Scalar>(
        Eigen::AngleAxis<Scalar>(rot_angle, rot_axis));

    std::string line;
    while (std::getline(ifs, line)) {
      std::vector<std::string> buf = stringutils::tokenize(line, ' ');

      if (buf.size() < 4 || buf[0] != "v") continue;

      Vec4x pos = Vec4x::Zero();
      for (int i = 0; i < 3; ++i)
        stringutils::extractFromString(buf[i + 1], pos(i));

      pos.segment<3>(0) = qrot * pos.segment<3>(0) + translate;

      m_particles.push_back(pos);
      m_fixed.push_back(0);
      m_surf_height.push_back(0.0);
      m_particle_groups.push_back(group);
      m_particle_closed_face_index.push_back(-1);
    }

    ifs.close();

    maxgroup = std::max(maxgroup, group);

    m_obj_reader_base.push_back(base_idx);
  }
}

void XMLReader::loadStrandParameters(rapidxml::xml_node<>* node,
                                     const double& dt) {
  rapidxml::xml_node<>* nd;

  const Scalar air_drag = m_simulation_params.m_airDrag;

  int paramsCount = 0;
  for (nd = node->first_node("StrandParameters"); nd;
       nd = nd->next_sibling("StrandParameters")) {
    // default values:
    double radius = 0.018;
    double YoungsModulus = 6.687e5;
    double shearModulus = 2.476e5;
    double density = 1.3;
    double viscosity = 0.0;
    double baseRotation = 0.;
    double stretchMultiplier = 1.0;
    double fixingMultiplier = 1.0;
    double minBendingAngle = 45.0 * M_PI / 180.0;
    double maxFlowGradientRatio = 10.0;
    double minBendingMultiplier = 1.0;
    double poissonRatio = 0.35;

    rapidxml::xml_node<>* subnd;

    loadParam(nd, "radius", radius);

    double biradius = radius;
    loadParam(nd, "biradius", biradius);
    loadParam(nd, "youngsModulus", YoungsModulus);
    loadParam(nd, "poissonRatio", poissonRatio);
    shearModulus = YoungsModulus / ((1.0 + poissonRatio) * 2.0);

    loadParam(nd, "density", density);
    loadParam(nd, "viscosity", viscosity);
    loadParam(nd, "baseRotation", baseRotation);
    loadParam(nd, "stretchMultiplier", stretchMultiplier);
    loadParam(nd, "fixingMultiplier", fixingMultiplier);
    loadParam(
        nd, "minBendingAngle", minBendingAngle,
        std::function<double(const double&)>(
            [](const double& val) -> double { return val * M_PI / 180.0; }));
    loadParam(nd, "minBendingMultiplier", minBendingMultiplier);
    loadParam(nd, "maxFlowGradientRatio", maxFlowGradientRatio);

    m_strand_parameters.push_back(ElasticStrandParameters(
        paramsCount, radius, biradius, YoungsModulus, shearModulus, density,
        viscosity, air_drag, baseRotation, stretchMultiplier, fixingMultiplier,
        minBendingAngle, maxFlowGradientRatio, minBendingMultiplier));

    ++paramsCount;
  }

  m_collision_parameters.resize(m_strand_parameters.size());

  int i = 0;
  for (CollisionParameters& param : m_collision_parameters) {
    setRodCollisionParameters(param, m_strand_parameters[i++]);
  }
}

template <typename T>
void XMLReader::loadParam(rapidxml::xml_node<>* nd, const char* name, T& param,
                          std::function<T(const T&)> const& post_process_func) {
  rapidxml::xml_node<>* subnd;
  if ((subnd = nd->first_node(name))) {
    std::string attribute(subnd->first_attribute("value")->value());
    if (!stringutils::extractFromString(attribute, param)) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of " << name << " attribute for "
                << nd->name() << ". Exiting." << std::endl;
      exit(1);
    }
    if (post_process_func) param = post_process_func(param);
  }
}

void XMLReader::setSimulationParameters() {
  m_simulation_params.m_statGathering = false;
#if defined(_OPENMP)
  m_simulation_params.m_numberOfThreads = omp_get_max_threads();
#else
  m_simulation_params.m_numberOfThreads = 1;
#endif
  m_simulation_params.m_rodSubSteps = 1;
  m_simulation_params.m_logLevel = MsgInfo::kInfo;
  m_simulation_params.m_useExactJacobian = false;
  m_simulation_params.m_useProjectedJacobian = false;
  m_simulation_params.m_useLengthProjection = true;
  m_simulation_params.m_usePreFilterGeometry = false;
  m_simulation_params.m_useApproxRodElasticFriction = true;
  m_simulation_params.m_skipRodMeshCollisions = false;
  m_simulation_params.m_skipRodRodCollisions = false;
  m_simulation_params.m_skipFlowFrictions = true;
  m_simulation_params.m_useCTRodRodCollisions = false;
  m_simulation_params.m_percentCTRodRodCollisionsAccept = 100.;
  m_simulation_params.m_useNonlinearContacts = false;
  m_simulation_params.m_solveLiquids = true;
  m_simulation_params.m_useAdditionalExternalFailSafe = false;
  m_simulation_params.m_useImpulseMethod = false;
  m_simulation_params.m_maxNewtonIterations = 10000;
  m_simulation_params.m_simulationManager_limitedMemory = false;
  m_simulation_params.m_gaussSeidelIterations = 75;
  m_simulation_params.m_gaussSeidelTolerance = 1e-4;
  m_simulation_params.m_pruneSelfCollisions = true;
  m_simulation_params.m_pruneExternalCollisions = true;
  m_simulation_params.m_stochasticPruning = 0.8;
  m_simulation_params.m_useDeterministicSolver = true;
  m_simulation_params.m_useSoftAttachConstraints = false;
  m_simulation_params.m_subSteps = 1;
  m_simulation_params.m_hairHairFrictionCoefficient = 0.3;
  m_simulation_params.m_hairMeshFrictionCoefficient = 0.0;
  m_simulation_params.m_airDrag = 0.0003;

  m_simulation_params.m_bogusAlgorithm =
      bogus::MecheFrictionProblem::ProjectedGradient;

  if (!m_scene_node) return;

  rapidxml::xml_node<>* nd = m_scene_node->first_node("SimulationParameters");

  if (!nd) return;

  loadParam(nd, "statGathering", m_simulation_params.m_statGathering);
  loadParam(nd, "numberOfThreads", m_simulation_params.m_numberOfThreads);
  loadParam(nd, "rodSubSteps", m_simulation_params.m_rodSubSteps);
  loadParam(nd, "logLevel", (int&)m_simulation_params.m_logLevel);
  loadParam(nd, "useLengthProjection",
            m_simulation_params.m_useLengthProjection);
  loadParam(nd, "usePreFilterGeometry",
            m_simulation_params.m_usePreFilterGeometry);
  loadParam(nd, "useExactJacobian", m_simulation_params.m_useExactJacobian);
  loadParam(nd, "useProjectedJacobian",
            m_simulation_params.m_useProjectedJacobian);
  loadParam(nd, "useApproxRodElasticFriction",
            m_simulation_params.m_useApproxRodElasticFriction);
  loadParam(nd, "skipRodMeshCollisions",
            m_simulation_params.m_skipRodMeshCollisions);
  loadParam(nd, "skipRodRodCollisions",
            m_simulation_params.m_skipRodRodCollisions);
  loadParam(nd, "skipFlowFrictions", m_simulation_params.m_skipFlowFrictions);
  loadParam(nd, "useCTRodRodCollisions",
            m_simulation_params.m_useCTRodRodCollisions);
  loadParam(nd, "percentCTRodRodCollisionsAccept",
            m_simulation_params.m_percentCTRodRodCollisionsAccept);
  loadParam(nd, "useNonlinearContacts",
            m_simulation_params.m_useNonlinearContacts);
  loadParam(nd, "solveLiquids", m_simulation_params.m_solveLiquids);
  loadParam(nd, "useAdditionalExternalFailSafe",
            m_simulation_params.m_useAdditionalExternalFailSafe);
  loadParam(nd, "useImpulseMethod", m_simulation_params.m_useImpulseMethod);
  loadParam(nd, "maxNewtonIterations",
            m_simulation_params.m_maxNewtonIterations);
  loadParam(nd, "simulationManagerLimitedMemory",
            m_simulation_params.m_simulationManager_limitedMemory);
  loadParam(nd, "gaussSeidelIterations",
            m_simulation_params.m_gaussSeidelIterations);
  loadParam(nd, "gaussSeidelTolerance",
            m_simulation_params.m_gaussSeidelTolerance);
  loadParam(nd, "pruneSelfCollisions",
            m_simulation_params.m_pruneSelfCollisions);
  loadParam(nd, "pruneExternalCollisions",
            m_simulation_params.m_pruneExternalCollisions);
  loadParam(nd, "stochasticPruning", m_simulation_params.m_stochasticPruning);
  loadParam(nd, "useDeterministicSolver",
            m_simulation_params.m_useDeterministicSolver);
  loadParam(nd, "useSoftAttachConstraints",
            m_simulation_params.m_useSoftAttachConstraints);
  loadParam(nd, "hairHairFrictionCoefficient",
            m_simulation_params.m_hairHairFrictionCoefficient);
  loadParam(nd, "hairMeshFrictionCoefficient",
            m_simulation_params.m_hairMeshFrictionCoefficient);
  loadParam(nd, "airDrag", m_simulation_params.m_airDrag);
  loadParam(nd, "subSteps", m_simulation_params.m_subSteps);

  rapidxml::xml_node<>* subnd;
  if ((subnd = nd->first_node("bogusAlgorithm"))) {
    std::string attribute(subnd->first_attribute("value")->value());
    if (attribute == "projectedgradient")
      m_simulation_params.m_bogusAlgorithm =
          bogus::MecheFrictionProblem::ProjectedGradient;
    else if (attribute == "matrixfreegaussseidel")
      m_simulation_params.m_bogusAlgorithm =
          bogus::MecheFrictionProblem::MatrixFreeGaussSeidel;
    else if (attribute == "gaussseidel")
      m_simulation_params.m_bogusAlgorithm =
          bogus::MecheFrictionProblem::GaussSeidel;
    else {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Invalid simulation parameter 'bogusAlgorithm' specified. "
                   "Exiting."
                << std::endl;
      exit(1);
    }
  }
}

void XMLReader::setRodCollisionParameters(
    CollisionParameters& param, const ElasticStrandParameters& strand_param) {
  param.m_externalCollisionsRadius = 0.005;
  param.m_selfCollisionsRadius = 0.005;
  param.m_constantCollisionRadius = true;
  param.m_maxNumCollisionsPerEdge = 6;
  param.m_impulseMaxNorm = 0.0;
  param.m_reactsToSelfCollisions = true;
  param.m_createsSelfCollisions = true;
  param.m_fakeLayering = true;
  param.m_frictionCoefficient =
      m_simulation_params.m_hairHairFrictionCoefficient;
  param.m_meshFrictionCoefficient =
      m_simulation_params.m_hairMeshFrictionCoefficient;

  if (!m_scene_node) return;

  rapidxml::xml_node<>* nd = m_scene_node->first_node("CollisionParameters");

  if (!nd) return;

  loadParam(nd, "externalCollisionsRadius", param.m_externalCollisionsRadius);
  loadParam(nd, "selfCollisionsRadius", param.m_selfCollisionsRadius);
  loadParam(nd, "constantCollisionRadius", param.m_constantCollisionRadius);
  loadParam(nd, "maxNumCollisionsPerEdge", param.m_maxNumCollisionsPerEdge);
  loadParam(nd, "impulseMaxNorm", param.m_impulseMaxNorm);
  loadParam(nd, "reactsToSelfCollisions", param.m_reactsToSelfCollisions);
  loadParam(nd, "createsSelfCollisions", param.m_createsSelfCollisions);
  loadParam(nd, "fakeLayering", param.m_fakeLayering);
  loadParam(nd, "frictionCoefficient", param.m_frictionCoefficient);
  loadParam(nd, "meshFrictionCoefficient", param.m_meshFrictionCoefficient);

  param.setAssociatedStrandParameters(strand_param);
}

void XMLReader::loadCheckpoint(rapidxml::xml_node<>* node, int& current_frame,
                               int& current_checkpoint) {
  for (rapidxml::xml_node<>* nd = node->first_node("checkpoint"); nd;
       nd = nd->next_sibling("checkpoint")) {
    DumpDataBinary data;
    if (nd->first_attribute("path")) {
      std::string file_path(nd->first_attribute("path")->value());

      ifstream ifs(file_path.c_str(), std::ios::binary);

      if (ifs.fail()) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to read file: " << file_path << ". Exiting."
                  << std::endl;
        exit(1);
      }
      // read general info
      ifs.read((char*)&(data.m_time), sizeof(Scalar));
      ifs.read((char*)&(data.current_frame), sizeof(int));
      ifs.read((char*)&(data.current_checkpoint), sizeof(int));
      ifs.read((char*)&(data.num_components), sizeof(int));

      // read fluid
      int num_particles;
      ifs.read((char*)&(num_particles), sizeof(int));

      if (num_particles) {
        data.m_x.resize(num_particles * 3);
        data.m_v.resize(num_particles * 3);
        data.m_m.resize(num_particles * 3);
        data.m_radius.resize(num_particles);
        data.m_J.resize(num_particles);
        data.m_vol.resize(num_particles);
        data.m_rest_vol.resize(num_particles);
        data.m_particle_group.resize(num_particles);
        data.m_classifier.resize(num_particles);
        data.m_weakened.resize(num_particles);
        data.m_components.resize(num_particles * data.num_components);
        data.m_proj_func.resize(num_particles * 3);

        data.m_Fe.resize(num_particles * 3, 3);
        data.m_b.resize(num_particles * 3, 3);
        data.m_B.resize(num_particles * 3, 3);
        data.m_b_trial.resize(num_particles * 3, 3);

        ifs.read((char*)data.m_x.data(), num_particles * 3 * sizeof(Scalar));
        ifs.read((char*)data.m_v.data(), num_particles * 3 * sizeof(Scalar));
        ifs.read((char*)data.m_m.data(), num_particles * 3 * sizeof(Scalar));
        ifs.read((char*)data.m_radius.data(), num_particles * sizeof(Scalar));
        ifs.read((char*)data.m_J.data(), num_particles * sizeof(Scalar));
        ifs.read((char*)data.m_vol.data(), num_particles * sizeof(Scalar));
        ifs.read((char*)data.m_rest_vol.data(), num_particles * sizeof(Scalar));
        ifs.read((char*)&(data.m_particle_group[0]),
                 num_particles * sizeof(int));
        ifs.read((char*)&(data.m_classifier[0]),
                 num_particles * sizeof(ParticleClassifier));
        ifs.read((char*)data.m_weakened.data(),
                 num_particles * sizeof(unsigned char));
        ifs.read((char*)data.m_components.data(),
                 num_particles * data.num_components * sizeof(Scalar));
        ifs.read((char*)data.m_proj_func.data(),
                 num_particles * 3 * sizeof(Scalar));

        ifs.read((char*)data.m_Fe.data(), num_particles * 9 * sizeof(Scalar));
        ifs.read((char*)data.m_b.data(), num_particles * 9 * sizeof(Scalar));
        ifs.read((char*)data.m_B.data(), num_particles * 9 * sizeof(Scalar));
        ifs.read((char*)data.m_b_trial.data(),
                 num_particles * 9 * sizeof(Scalar));
      }

      // read hairs
      int num_strands;
      ifs.read((char*)&(num_strands), sizeof(int));

      data.m_currentDOFs.resize(num_strands);
      data.m_currentAreaDOFs.resize(num_strands);
      data.m_velocities.resize(num_strands);
      data.m_flow_velocities.resize(num_strands);
      data.m_flow_strain.resize(num_strands);
      data.m_flow_components.resize(num_strands);
      data.m_flow_reservoirs.resize(num_strands);

      for (int i = 0; i < num_strands; ++i) {
        int num_verts;
        ifs.read((char*)&(num_verts), sizeof(int));

        data.m_currentDOFs[i].resize(num_verts * 4 - 1);
        data.m_currentAreaDOFs[i].resize(num_verts);
        data.m_velocities[i].resize(num_verts * 4 - 1);
        data.m_flow_velocities[i].resize(num_verts - 1);
        data.m_flow_strain[i].resize(num_verts);
        data.m_flow_components[i].resize(num_verts * data.num_components);

        ifs.read((char*)(data.m_currentDOFs[i].data()),
                 data.m_currentDOFs[i].size() * sizeof(Scalar));
        ifs.read((char*)(data.m_currentAreaDOFs[i].data()),
                 data.m_currentAreaDOFs[i].size() * sizeof(Scalar));
        ifs.read((char*)(data.m_velocities[i].data()),
                 data.m_velocities[i].size() * sizeof(Scalar));
        ifs.read((char*)(data.m_flow_velocities[i].data()),
                 data.m_flow_velocities[i].size() * sizeof(Scalar));
        ifs.read((char*)(data.m_flow_strain[i].data()),
                 data.m_flow_strain[i].size() * sizeof(Scalar));
        ifs.read((char*)(data.m_flow_components[i].data()),
                 data.m_flow_components[i].size() * sizeof(Scalar));
        ifs.read((char*)(data.m_flow_reservoirs[i].data()), 2 * sizeof(Scalar));
      }

      // read goals
      int num_goals;
      ifs.read((char*)&(num_goals), sizeof(int));
      data.m_strand_goals.resize(num_goals);
      if (num_goals) {
        ifs.read((char*)&(data.m_strand_goals[0]),
                 sizeof(std::pair<Vec2i, Vec4x>) * num_goals);
      }

      // read groups
      int num_groups;
      ifs.read((char*)&num_groups, sizeof(int));
      if (num_groups) {
        data.m_group_pos.resize(num_groups);
        data.m_group_scale.resize(num_groups);
        data.m_group_rot.resize(num_groups);

        ifs.read((char*)&(data.m_group_pos[0]), sizeof(Vec3x) * num_groups);
        ifs.read((char*)&(data.m_group_scale[0]), sizeof(Vec3x) * num_groups);
        ifs.read((char*)&(data.m_group_rot[0]),
                 sizeof(Eigen::Quaternion<Scalar>) * num_groups);
      }

      // read fields
      int num_fields;
      ifs.read((char*)&num_fields, sizeof(int));

      if (num_fields) {
        data.m_field_i_frames.resize(num_fields);
        data.m_field_center.resize(num_fields);
        data.m_field_rot.resize(num_fields);
        data.m_field_future_scale.resize(num_fields);

        ifs.read((char*)&(data.m_field_i_frames[0]), sizeof(int) * num_fields);
        ifs.read((char*)&(data.m_field_center[0]), sizeof(Vec3x) * num_fields);
        ifs.read((char*)&(data.m_field_rot[0]),
                 sizeof(Eigen::Quaternion<Scalar>) * num_fields);
        ifs.read((char*)&(data.m_field_future_scale[0]),
                 sizeof(Vec3x) * num_fields);
      }

      data.m_field_parameter.resize(num_fields);

      for (int i = 0; i < num_fields; ++i) {
        int num_params;
        ifs.read((char*)&(num_params), sizeof(int));

        data.m_field_parameter[i].resize(num_params);
        ifs.read((char*)(data.m_field_parameter[i].data()),
                 num_params * sizeof(Scalar));
      }

      // write mesh
      int num_meshes;
      ifs.read((char*)&num_meshes, sizeof(int));
      data.m_mesh_current_frames.resize(num_meshes);
      data.m_previous_mesh_vertices.resize(num_meshes);
      data.m_current_mesh_vertices.resize(num_meshes);
      data.m_next_mesh_vertices.resize(num_meshes);
      data.m_current_mesh_displacement.resize(num_meshes);

      for (int i = 0; i < num_meshes; ++i) {
        int num_verts;
        ifs.read((char*)&num_verts, sizeof(int));
        data.m_previous_mesh_vertices[i].resize(num_verts);
        data.m_current_mesh_vertices[i].resize(num_verts);
        data.m_next_mesh_vertices[i].resize(num_verts);
        data.m_current_mesh_displacement[i].resize(num_verts);

        ifs.read((char*)&(data.m_mesh_current_frames[i]), sizeof(int));
        ifs.read((char*)&(data.m_previous_mesh_vertices[i][0]),
                 sizeof(Vec3x) * data.m_previous_mesh_vertices[i].size());
        ifs.read((char*)&(data.m_current_mesh_vertices[i][0]),
                 sizeof(Vec3x) * data.m_current_mesh_vertices[i].size());
        ifs.read((char*)&(data.m_next_mesh_vertices[i][0]),
                 sizeof(Vec3x) * data.m_next_mesh_vertices[i].size());
        ifs.read((char*)&(data.m_current_mesh_displacement[i][0]),
                 sizeof(Vec3x) * data.m_current_mesh_displacement[i].size());
      }

      ifs.close();

      // set general
      m_t = data.m_time;
      m_stepper->setTime(data.m_time);
      current_frame = data.current_frame;
      current_checkpoint = data.current_checkpoint;

      // set fluid
      if (m_fluidScripting_controllers.size() &&
          m_fluidScripting_controllers[0]) {
        std::shared_ptr<FluidScriptingController> controller =
            m_fluidScripting_controllers[0];

        controller->conservativeResizeParticles(num_particles);
        controller->getX() = data.m_x;
        controller->getV() = data.m_v;
        controller->getM() = data.m_m;
        controller->getRadius() = data.m_radius;
        controller->getJ() = data.m_J;
        controller->getVol() = data.m_vol;
        controller->getRestVol() = data.m_rest_vol;
        controller->getParticleGroup() = data.m_particle_group;
        controller->getClassifier() = data.m_classifier;
        controller->getWeakened() = data.m_weakened;
        controller->getComponents() = data.m_components;
        controller->getProjFunc() = data.m_proj_func;

        controller->getFe() = data.m_Fe;
        controller->getFePlus() = data.m_Fe;
        controller->getb() = data.m_b;
        controller->getB() = data.m_B;
        controller->getbtrial() = data.m_b_trial;
        controller->getbPlus() = data.m_b;
      }

      // set hairs
      assert((int)m_strands.size() == num_strands);
      for (int i = 0; i < num_strands; ++i) {
        m_strands[i]->setCurrentDegreesOfFreedom(data.m_currentDOFs[i]);
        m_strands[i]->setFutureDegreesOfFreedom(data.m_currentDOFs[i] -
                                                data.m_velocities[i] * m_dt);
        m_strands[i]->setCurrentAreaDegreesOfFreedom(data.m_currentAreaDOFs[i]);
        m_strands[i]->setFutureAreaDegreesOfFreedom(data.m_currentAreaDOFs[i]);
        m_strands[i]->getStepper()->velocities() = data.m_velocities[i];
        m_strands[i]->getStepper()->newVelocities() = data.m_velocities[i];
        m_strands[i]->dynamics().getDisplacements() =
            data.m_velocities[i] * m_dt;
        m_strands[i]->getStepper()->flowVelocities() =
            data.m_flow_velocities[i];
        m_strands[i]->getStepper()->newFlowVelocities() =
            data.m_flow_velocities[i];
        m_strands[i]->getStepper()->flowStrain() = data.m_flow_strain[i];
        m_strands[i]->getStepper()->flowNewStrain() = data.m_flow_strain[i];
        m_strands[i]->getStepper()->flowComponents() =
            data.m_flow_components[i];
        m_strands[i]->getStepper()->flowNewComponents() =
            data.m_flow_components[i];
        m_strands[i]->getReservoir() = data.m_flow_reservoirs[i];

        const int num_verts = m_strands[i]->getNumVertices();
        for (int k = 0; k < num_verts; ++k) {
          if (m_strands[i]->isVertexFreezed(k)) {
            m_strands[i]->dynamics().getScriptingController()->freezeVertices(
                k);
          }

          if (k != num_verts - 1 && m_strands[i]->isThetaFreezed(k)) {
            m_strands[i]->dynamics().getScriptingController()->freezeTheta(k);
          }
        }
      }

      for (const std::pair<Vec2i, Vec4x>& goal : data.m_strand_goals) {
        if (m_strands[goal.first[0]]->isVertexGoaled(goal.first[1]))
          m_strands[goal.first[0]]
              ->dynamics()
              .getScriptingController()
              ->setVertexGoal(
                  goal.first[1],
                  Vec3x(goal.second(0), goal.second(1), goal.second(2)));

        if (m_strands[goal.first[0]]->isThetaGoaled(goal.first[1]))
          m_strands[goal.first[0]]
              ->dynamics()
              .getScriptingController()
              ->setThetaGoal(goal.first[1], goal.second(3));
      }

      // set groups
      assert(m_groups.size() >= data.m_group_pos.size());

      for (int i = 0;
           i < (int)m_group_pos.size() && i < (int)data.m_group_pos.size();
           ++i) {
        m_group_prev_pos[i] = m_group_pos[i] = data.m_group_pos[i];
        m_group_prev_scale[i] = m_group_scale[i] = data.m_group_scale[i];
        m_group_prev_rot[i] = m_group_rot[i] = data.m_group_rot[i];
      }

      // set fields
      assert(m_fields.size() + m_sources.size() + m_terminators.size() >=
             data.m_field_i_frames.size());

      int k = 0;
      for (int i = 0;
           i < (int)m_fields.size() && k < (int)data.m_field_i_frames.size();
           ++i, ++k) {
        m_fields[i].i_frame = data.m_field_i_frames[k];
        m_fields[i].future_center = m_fields[i].center = data.m_field_center[k];
        m_fields[i].future_rot = m_fields[i].rot = data.m_field_rot[k];
        m_fields[i].future_scale = data.m_field_future_scale[k];
        m_fields[i].parameter = data.m_field_parameter[k];
      }

      for (int i = 0;
           i < (int)m_sources.size() && k < (int)data.m_field_i_frames.size();
           ++i, ++k) {
        m_sources[i].i_frame = data.m_field_i_frames[k];
        m_sources[i].future_center = m_sources[i].center =
            data.m_field_center[k];
        m_sources[i].future_rot = m_sources[i].rot = data.m_field_rot[k];
        m_sources[i].future_scale = data.m_field_future_scale[k];
        m_sources[i].parameter = data.m_field_parameter[k];
      }

      for (int i = 0; i < (int)m_terminators.size() &&
                      k < (int)data.m_field_i_frames.size();
           ++i, ++k) {
        m_terminators[i].i_frame = data.m_field_i_frames[k];
        m_terminators[i].future_center = m_terminators[i].center =
            data.m_field_center[k];
        m_terminators[i].future_rot = m_terminators[i].rot =
            data.m_field_rot[k];
        m_terminators[i].future_scale = data.m_field_future_scale[k];
        m_terminators[i].parameter = data.m_field_parameter[k];
      }

      // set mesh
      assert(m_meshScripting_controllers.size() == num_meshes);

      for (int i = 0; i < num_meshes; ++i) {
        m_meshScripting_controllers[i]->setIFrame(
            data.m_mesh_current_frames[i]);
        m_meshScripting_controllers[i]->getPreviousMesh()->m_vertices =
            data.m_previous_mesh_vertices[i];
        m_meshScripting_controllers[i]->getCurrentMesh()->m_vertices =
            data.m_current_mesh_vertices[i];
        m_meshScripting_controllers[i]->getNextMesh()->m_vertices =
            data.m_next_mesh_vertices[i];
        m_meshScripting_controllers[i]
            ->getCurrentMesh()
            ->m_previous_displacements =
            m_meshScripting_controllers[i]->getCurrentMesh()->m_displacements =
                data.m_current_mesh_displacement[i];

        m_meshScripting_controllers[i]->updateMeshNormalArea();
      }
    }
  }

  executeCallback();
}

void XMLReader::loadHairPose(rapidxml::xml_node<>* node) {
  for (rapidxml::xml_node<>* nd = node->first_node("pose"); nd;
       nd = nd->next_sibling("pose")) {
    HairPose pose;
    pose.loaded = false;
    pose.time = 0.0;

    if (nd->first_attribute("time")) {
      std::string attribute(nd->first_attribute("time")->value());
      if (!stringutils::extractFromString(attribute, pose.time)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of time attribute for pose. Value "
                     "must be double. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    std::unordered_map<int, int> strand_indices_selection;

    for (rapidxml::xml_node<>* subnd = nd->first_node("selection"); subnd;
         subnd = subnd->next_sibling("selection")) {
      if (subnd->first_attribute("index")) {
        std::string attribute(subnd->first_attribute("index")->value());

        std::vector<string> data;
        stringutils::split(attribute, ' ', data);

        for (const string& s : data) {
          int si;
          if (stringutils::extractFromString(s, si)) {
            int mapped_i = strand_indices_selection.size();
            strand_indices_selection[si] = mapped_i;
          }
        }
      }
    }

    Vec3x center = Vec3x::Zero();
    if (nd->first_attribute("cx")) {
      std::string attribute(nd->first_attribute("cx")->value());
      if (!stringutils::extractFromString(attribute, center(0))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of center(0) attribute for pose "
                     "parameters. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("cy")) {
      std::string attribute(nd->first_attribute("cy")->value());
      if (!stringutils::extractFromString(attribute, center(1))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of center(1) attribute for pose "
                     "parameters. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("cz")) {
      std::string attribute(nd->first_attribute("cz")->value());
      if (!stringutils::extractFromString(attribute, center(2))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of center(2) attribute for pose "
                     "parameters. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("path")) {
      std::string file_path(nd->first_attribute("path")->value());

      ifstream ifs(file_path.c_str());

      if (ifs.fail()) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to read file: " << file_path << ". Exiting."
                  << std::endl;
        exit(1);
      }

      string line;
      bool start_data = false;
      std::vector<Vec4x> pos;
      std::vector<Scalar> areas;
      std::vector<std::pair<int, int> > segment;
      std::unordered_map<int, int> seg_max_inds;

      int num_segs = 0;

      int idx_x = 0;
      int idx_y = 1;
      int idx_z = 2;
      int idx_theta = 3;
      int idx_seg = 4;
      int idx_ra = 5;
      int idx_ha = 7;
      int idx_actual = 10;

      int prop_order = 0;

      while (getline(ifs, line)) {
        vector<string> data;
        stringutils::split(line, ' ', data);

        if (start_data) {
          if (data.size() < 8) continue;

          Vec4x p;
          Scalar r, h, area;
          int seg;
          int actual;

          stringutils::extractFromString(data[idx_actual], actual);
          if (!actual) continue;

          stringutils::extractFromString(data[idx_x], p(0));
          stringutils::extractFromString(data[idx_y], p(1));
          stringutils::extractFromString(data[idx_z], p(2));
          stringutils::extractFromString(data[idx_theta], p(3));
          stringutils::extractFromString(data[idx_seg], seg);
          stringutils::extractFromString(data[idx_ra], r);
          stringutils::extractFromString(data[idx_ha], h);

          if (strand_indices_selection.size()) {
            if (strand_indices_selection.find(seg) ==
                strand_indices_selection.end())
              continue;
            else
              seg = strand_indices_selection[seg];
          }

          h -= r;
          area = M_PI * (h + 2.0 * r) * h;

          p.segment<3>(0) += center;

          pos.push_back(p);
          areas.push_back(area);

          num_segs = std::max(num_segs, seg + 1);
          if (seg_max_inds.find(seg) == seg_max_inds.end()) {
            seg_max_inds[seg] = 1;
          } else {
            seg_max_inds[seg]++;
          }

          segment.push_back(std::pair<int, int>(seg, seg_max_inds[seg] - 1));

        } else if (line.substr(0, 10) == "end_header") {
          start_data = true;
        } else {
          if (data.size() == 0) continue;

          if (data[0] == "property") {
            if (data.size() != 3) continue;

            if (data[2] == "x")
              idx_x = prop_order;
            else if (data[2] == "y")
              idx_y = prop_order;
            else if (data[2] == "z")
              idx_z = prop_order;
            else if (data[2] == "theta")
              idx_theta = prop_order;
            else if (data[2] == "segment")
              idx_seg = prop_order;
            else if (data[2] == "ra")
              idx_ra = prop_order;
            else if (data[2] == "ha")
              idx_ha = prop_order;
            else if (data[2] == "actual")
              idx_actual = prop_order;

            prop_order++;
          } else if (data[0] == "comment") {
            continue;
          }
        }
      }

      pose.dofs.resize(num_segs);
      pose.areas.resize(num_segs);

      for (int i = 0; i < num_segs; ++i) {
        const int num_verts = seg_max_inds[i];
        pose.dofs[i].resize(num_verts * 4 - 1);
        pose.areas[i].resize(num_verts);
      }

      const int total_num_verts = pos.size();
      for (int i = 0; i < total_num_verts; ++i) {
        const int num_vert_strand = seg_max_inds[segment[i].first];
        if (segment[i].second == num_vert_strand - 1) {
          pose.dofs[segment[i].first].segment<3>(segment[i].second * 4) =
              pos[i].segment<3>(0);
        } else {
          pose.dofs[segment[i].first].segment<4>(segment[i].second * 4) =
              pos[i];
        }

        pose.areas[segment[i].first][segment[i].second] = areas[i];
      }

      m_hairposes.emplace_back(pose);

      ifs.close();
    } else {
      continue;
    }
  }

  executeHairPose();
  executeCallback();
}

void XMLReader::loadScripts(rapidxml::xml_node<>* node) {
  for (rapidxml::xml_node<>* nd = node->first_node("script"); nd;
       nd = nd->next_sibling("script")) {
    Script scr;
    scr.used = false;

    rapidxml::xml_attribute<>* typend = NULL;
    typend = nd->first_attribute("type");
    if (typend) {
      std::string handlertype(typend->value());
      if (handlertype == "rotate")
        scr.type = Script::ROTATE;
      else if (handlertype == "translate")
        scr.type = Script::TRANSLATE;
      else if (handlertype == "swirl")
        scr.type = Script::SWIRL;
      else if (handlertype == "scale")
        scr.type = Script::SCALE;
      else if (handlertype == "switch")
        scr.type = Script::SWITCH;
      else {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Invalid script 'type' attribute specified. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Invalid script 'type' attribute specified. Exiting."
                << std::endl;
      exit(1);
    }

    scr.func = Script::CUBIC;
    typend = nd->first_attribute("func");
    if (typend) {
      std::string handlertype(typend->value());
      if (handlertype == "cubic")
        scr.func = Script::CUBIC;
      else if (handlertype == "cosine")
        scr.func = Script::COSINE;
      else if (handlertype == "weno")
        scr.func = Script::WENO;
      else {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Invalid script 'func' attribute specified. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    scr.base_pos = 0.0;

    if (nd->first_attribute("x")) {
      std::string attribute(nd->first_attribute("x")->value());
      if (!stringutils::extractFromString(attribute, scr.v(0))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of x attribute for script. Value "
                     "must be double. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("y")) {
      std::string attribute(nd->first_attribute("y")->value());
      if (!stringutils::extractFromString(attribute, scr.v(1))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of y attribute for script. Value "
                     "must be double. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("z")) {
      std::string attribute(nd->first_attribute("z")->value());
      if (!stringutils::extractFromString(attribute, scr.v(2))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of z attribute for script. Value "
                     "must be double. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("w")) {
      std::string attribute(nd->first_attribute("w")->value());
      if (!stringutils::extractFromString(attribute, scr.v(3))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of w attribute for script. Value "
                     "must be double. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.transform_with_origin = false;
    if (nd->first_attribute("ox")) {
      std::string attribute(nd->first_attribute("ox")->value());
      if (!stringutils::extractFromString(attribute, scr.origin(0))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of ox attribute for script. Value "
                     "must be double. Exiting."
                  << std::endl;
        exit(1);
      }
      scr.transform_with_origin = true;
    }
    if (nd->first_attribute("oy")) {
      std::string attribute(nd->first_attribute("oy")->value());
      if (!stringutils::extractFromString(attribute, scr.origin(1))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of oy attribute for script. Value "
                     "must be double. Exiting."
                  << std::endl;
        exit(1);
      }
      scr.transform_with_origin = true;
    }
    if (nd->first_attribute("oz")) {
      std::string attribute(nd->first_attribute("oz")->value());
      if (!stringutils::extractFromString(attribute, scr.origin(2))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of oz attribute for script. Value "
                     "must be double. Exiting."
                  << std::endl;
        exit(1);
      }
      scr.transform_with_origin = true;
    }

    scr.switch_mask = 0;
    if (nd->first_attribute("mask")) {
      std::string attribute(nd->first_attribute("mask")->value());
      if (!stringutils::extractFromString(attribute, scr.switch_mask)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of mask attribute for script. "
                     "Value must be int. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.start = 0.0;
    if (nd->first_attribute("start")) {
      std::string attribute(nd->first_attribute("start")->value());
      if (!stringutils::extractFromString(attribute, scr.start)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of start attribute for script. "
                     "Value must be double. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.end = m_end_time;
    if (nd->first_attribute("end")) {
      std::string attribute(nd->first_attribute("end")->value());
      if (!stringutils::extractFromString(attribute, scr.end)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of end attribute for script. "
                     "Value must be double. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.ease_start = scr.ease_end = (scr.end - scr.start) / 3.0;

    if (nd->first_attribute("easestart")) {
      std::string attribute(nd->first_attribute("easestart")->value());
      if (!stringutils::extractFromString(attribute, scr.ease_start)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of start attribute for script. "
                     "Value must be double. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("easeend")) {
      std::string attribute(nd->first_attribute("easeend")->value());
      if (!stringutils::extractFromString(attribute, scr.ease_end)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of start attribute for script. "
                     "Value must be double. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.amplitude = 1.0;
    if (nd->first_attribute("amplitude")) {
      std::string attribute(nd->first_attribute("amplitude")->value());
      if (!stringutils::extractFromString(attribute, scr.amplitude)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of amplitude attribute for "
                     "script. Value must be double. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.base_dt = 1.0;
    if (nd->first_attribute("dt")) {
      std::string attribute(nd->first_attribute("dt")->value());
      if (!stringutils::extractFromString(attribute, scr.base_dt)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of dt attribute for script. Value "
                     "must be double. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("base")) {
      std::string attribute(nd->first_attribute("base")->value());
      std::vector<std::string> bases = stringutils::split(attribute, ' ');

      scr.base_vertices.reserve(bases.size());
      for (const std::string& str : bases) {
        double y = 0;
        stringutils::extractFromString(str, y);
        scr.base_vertices.push_back(y);
      }
    }

    scr.frequency = 1.0;
    if (nd->first_attribute("frequency")) {
      std::string attribute(nd->first_attribute("frequency")->value());
      if (!stringutils::extractFromString(attribute, scr.frequency)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of frequency attribute for "
                     "script. Value must be double. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.group_index = 0;
    if (nd->first_attribute("group")) {
      std::string attribute(nd->first_attribute("group")->value());
      if (!stringutils::extractFromString(attribute, scr.group_index)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of i attribute for script. Value "
                     "must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.transform_global = false;
    if (nd->first_attribute("global")) {
      std::string attribute(nd->first_attribute("global")->value());
      if (!stringutils::extractFromString(attribute, scr.transform_global)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of global attribute for script. "
                     "Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.m_scene = this;
    m_scripts.push_back(scr);
  }
}

void XMLReader::loadIntegrator(rapidxml::xml_node<>* node, double& dt) {
  {
    rapidxml::xml_node<>* nd = node->first_node("integrator");
    if (nd == NULL) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " No integrator specified. Exiting." << std::endl;
      exit(1);
    }

    rapidxml::xml_attribute<>* dtnd = nd->first_attribute("dt");
    if (dtnd == NULL) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " No integrator 'dt' attribute specified. Exiting."
                << std::endl;
      exit(1);
    } else {
      std::string attribute(dtnd->value());
      if (!stringutils::extractFromString(attribute, dt)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of dt attribute for integrator. "
                     "Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("cfl")) {
      std::string attribute(nd->first_attribute("cfl")->value());
      if (!stringutils::extractFromString(attribute, m_cfl_number)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of cfl attribute for integrator. "
                     "Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }
  }

  {
    rapidxml::xml_node<>* nd = node->first_node("duration");
    if (nd == NULL) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " No duration specified. Exiting." << std::endl;
      exit(1);
    }

    rapidxml::xml_attribute<>* dtnd = nd->first_attribute("time");
    if (dtnd == NULL) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " No duration 'time' attribute specified. Exiting."
                << std::endl;
      exit(1);
    }

    m_end_time = 0.0;
    if (!stringutils::extractFromString(std::string(dtnd->value()),
                                        m_end_time)) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse 'time' attribute for duration. Value must "
                   "be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }
}

void XMLReader::loadHairs(rapidxml::xml_node<>* node) {
  assert(node != NULL);

  int num_components = 1;

  rapidxml::xml_node<>* nd = node->first_node("liquidinfo");
  if (nd) {
    loadParam(nd, "numComponents", num_components);
  }

  // Count the number of particles, edges, and strands
  int numstrands = 0;
  int numparticles = m_particles.size();

  m_global_to_local.assign(numparticles, std::pair<int, int>(-1, -1));
  m_local_to_global.resize(0);

  for (rapidxml::xml_node<>* nd = node->first_node("hair"); nd;
       nd = nd->next_sibling("hair")) {
    int paramsIndex = -1;
    if (nd->first_attribute("params")) {
      std::string attribute(nd->first_attribute("params")->value());
      if (!stringutils::extractFromString(attribute, paramsIndex)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of params (StrandParameters "
                     "index) for hair "
                  << numstrands << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of params (StrandParameters index) "
                   "for hair "
                << numstrands << ". Exiting." << std::endl;
      exit(1);
    }

    if (paramsIndex == -1) continue;

    int start = 0;
    if (nd->first_attribute("start")) {
      std::string attribute(nd->first_attribute("start")->value());
      if (!stringutils::extractFromString(attribute, start)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of start attribute for hair "
                  << numstrands << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    int count = 0;
    if (nd->first_attribute("count")) {
      std::string attribute(nd->first_attribute("count")->value());
      if (!stringutils::extractFromString(attribute, count)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of count attribute for hair "
                  << numstrands << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    std::vector<int> particle_indices;

    if (count == 0) {
      for (rapidxml::xml_node<>* subnd = nd->first_node("p"); subnd;
           subnd = subnd->next_sibling("p")) {
        int id = -1;
        if (subnd->first_attribute("i")) {
          std::string attribute(subnd->first_attribute("i")->value());
          if (!stringutils::extractFromString(attribute, id)) {
            std::cerr << outputmod::startred
                      << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                      << " Failed to parse value of count attribute for p id "
                      << numstrands << ". Value must be integer. Exiting."
                      << std::endl;
            exit(1);
          }
        } else {
          std::cerr << outputmod::startred
                    << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                    << " Failed to parse value of count attribute for p id "
                    << numstrands << ". No attribute id." << std::endl;
          exit(1);
        }
        particle_indices.push_back(id);
      }

      count = (int)particle_indices.size();
    } else {
      particle_indices.resize(count);

      for (int i = 0; i < count; ++i) {
        const int pidx = start + i;

        particle_indices[i] = pidx;
      }
    }

    int num_DoFs = 4 * count - 1;
    VecXx dofs(num_DoFs);
    dofs.setZero();

    VecXx area_dofs(count);
    area_dofs.setZero();

    DOFScriptingController* controller = new DOFScriptingController();

    for (int i = 0; i < count; ++i) {
      if (i == count - 1) {
        dofs.segment<3>(i * 4) = Vec3x(m_particles[particle_indices[i]](0),
                                       m_particles[particle_indices[i]](1),
                                       m_particles[particle_indices[i]](2));
      } else {
        dofs.segment<4>(i * 4) = m_particles[particle_indices[i]];
      }

      const Scalar ra = m_strand_parameters[paramsIndex].getRadiusA(i, count);
      const Scalar rb = m_strand_parameters[paramsIndex].getRadiusB(i, count);
      const Scalar h = m_surf_height[particle_indices[i]];

      area_dofs(i) = M_PI * h * (h + ra + rb);

      if (m_fixed[particle_indices[i]] & 1) {
        if (m_simulation_params.m_useSoftAttachConstraints) {
          Vec3x pos = m_particles[particle_indices[i]].segment<3>(0);
          controller->setVertexGoal(i, pos);
        } else {
          controller->freezeVertices(i);
        }
      }

      if (m_fixed[particle_indices[i]] & 2) {
        //                controller->freezeTheta(i);
        if (m_simulation_params.m_useSoftAttachConstraints) {
          controller->setThetaGoal(i, m_particles[particle_indices[i]](3));
        } else {
          controller->freezeTheta(i);
        }
      }

      if (m_fixed[particle_indices[i]]) {
        const int group_idx = m_particle_groups[particle_indices[i]];
        assert(group_idx >= 0 && group_idx < (int)m_groups.size());

        m_groups[group_idx].push_back(std::pair<int, int>(numstrands, i));
      }

      m_global_to_local[particle_indices[i]] =
          std::pair<int, int>(numstrands, i);
    }

    ElasticStrand* strand =
        new ElasticStrand(dofs, area_dofs, m_strand_parameters[paramsIndex],
                          m_collision_parameters[paramsIndex], controller);

    strand->setGlobalIndex(numstrands);
    m_strands.push_back(strand);

    RodData* rd = new RodData(*strand, *controller, num_components);
    m_rodDatum.push_back(rd);

    m_local_to_global.push_back(particle_indices);

    ++numstrands;
  }

  int i_obj_strand = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("hairobj"); nd;
       nd = nd->next_sibling("hairobj")) {
    std::string szfn;
    if (nd->first_attribute("filename")) {
      szfn = nd->first_attribute("filename")->value();
    }

    if (szfn.empty()) continue;

    int base_idx = m_obj_reader_base[i_obj_strand];

    int paramsIndex = 0;
    if (nd->first_attribute("params")) {
      std::string attribute(nd->first_attribute("params")->value());
      if (!stringutils::extractFromString(attribute, paramsIndex)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of params attribute for obj "
                  << i_obj_strand << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    struct setFlow {
      int start;
      int end;
      double h;
    };

    struct setFixed {
      int start;
      int end;
      unsigned fixed;
    };

    std::vector<setFlow> flow_settings;

    for (rapidxml::xml_node<>* subnd = nd->first_node("flow"); subnd;
         subnd = subnd->next_sibling("flow")) {
      setFlow s = {0, 0, 0.0};

      if (subnd->first_attribute("start")) {
        std::string attribute(subnd->first_attribute("start")->value());
        if (!stringutils::extractFromString(attribute, s.start)) {
          exit(1);
        }
      }

      if (subnd->first_attribute("end")) {
        std::string attribute(subnd->first_attribute("end")->value());
        if (!stringutils::extractFromString(attribute, s.end)) {
          exit(1);
        }
      }

      if (subnd->first_attribute("value")) {
        std::string attribute(subnd->first_attribute("value")->value());
        if (!stringutils::extractFromString(attribute, s.h)) {
          exit(1);
        }
      }

      flow_settings.push_back(s);
    }

    std::vector<setFixed> fix_settings;
    for (rapidxml::xml_node<>* subnd = nd->first_node("fixed"); subnd;
         subnd = subnd->next_sibling("fixed")) {
      setFixed s = {0, 0, 0U};

      if (subnd->first_attribute("start")) {
        std::string attribute(subnd->first_attribute("start")->value());
        if (!stringutils::extractFromString(attribute, s.start)) {
          exit(1);
        }
      }

      if (subnd->first_attribute("end")) {
        std::string attribute(subnd->first_attribute("end")->value());
        if (!stringutils::extractFromString(attribute, s.end)) {
          exit(1);
        }
      }

      if (subnd->first_attribute("value")) {
        std::string attribute(subnd->first_attribute("value")->value());
        if (!stringutils::extractFromString(attribute, s.fixed)) {
          exit(1);
        }
      }

      fix_settings.push_back(s);
    }

    std::unordered_set<int> strand_indices_selection;

    for (rapidxml::xml_node<>* subnd = nd->first_node("selection"); subnd;
         subnd = subnd->next_sibling("selection")) {
      if (subnd->first_attribute("index")) {
        std::string attribute(subnd->first_attribute("index")->value());

        std::vector<string> data;
        stringutils::split(attribute, ' ', data);

        for (const string& s : data) {
          int si;
          if (stringutils::extractFromString(s, si)) {
            strand_indices_selection.insert(si);
          }
        }
      }
    }

    std::string line;
    std::ifstream ifs(szfn);

    if (ifs.fail()) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to read file: " << szfn << ". Exiting."
                << std::endl;
      exit(1);
    }

    int linenumber = 0;

    while (getline(ifs, line)) {
      std::vector<std::string> tokens = stringutils::tokenize(line, ' ');
      if (tokens.size() < 1 || tokens[0] != "l") continue;

      if (strand_indices_selection.size() &&
          strand_indices_selection.find(linenumber) ==
              strand_indices_selection.end()) {
        linenumber++;
        continue;
      }

      std::vector<int> particle_indices;

      for (int i = 1; i < (int)tokens.size(); ++i) {
        int pidx = -1;
        stringutils::extractFromString(tokens[i], pidx);
        if (pidx == -1) continue;

        particle_indices.push_back(base_idx + pidx - 1);
      }

      std::sort(particle_indices.begin(), particle_indices.end());

      int count = (int)particle_indices.size();

      int num_DoFs = 4 * count - 1;
      VecXx dofs(num_DoFs);
      dofs.setZero();

      VecXx area_dofs(count);
      area_dofs.setZero();

      DOFScriptingController* controller = new DOFScriptingController();

      for (const setFlow& s : flow_settings) {
        for (int i = s.start; i <= s.end; ++i) {
          const Scalar ra =
              m_strand_parameters[paramsIndex].getRadiusA(i, count);
          const Scalar rb =
              m_strand_parameters[paramsIndex].getRadiusB(i, count);

          area_dofs(i) = M_PI * s.h * (s.h + ra + rb);
          m_surf_height[particle_indices[i]] = s.h;
        }
      }

      for (const setFixed& s : fix_settings) {
        for (int i = s.start; i <= s.end; ++i) {
          m_fixed[particle_indices[i]] = s.fixed;

          if (m_fixed[particle_indices[i]] & 1) {
            //                controller->freezeVertices(i);
            if (m_simulation_params.m_useSoftAttachConstraints) {
              Vec3x pos = m_particles[particle_indices[i]].segment<3>(0);
              controller->setVertexGoal(i, pos);
            } else {
              controller->freezeVertices(i);
            }
          }

          if (m_fixed[particle_indices[i]] & 2) {
            //                controller->freezeTheta(i);
            if (m_simulation_params.m_useSoftAttachConstraints) {
              controller->setThetaGoal(i, m_particles[particle_indices[i]](3));
            } else {
              controller->freezeTheta(i);
            }
          }

          if (m_fixed[particle_indices[i]]) {
            const int group_idx = m_particle_groups[particle_indices[i]];
            assert(group_idx >= 0 && group_idx < (int)m_groups.size());

            m_groups[group_idx].push_back(std::pair<int, int>(numstrands, i));
          }
        }
      }

      for (int i = 0; i < count; ++i) {
        m_global_to_local[particle_indices[i]] =
            std::pair<int, int>(numstrands, i);
        if (i == count - 1) {
          dofs.segment<3>(i * 4) = Vec3x(m_particles[particle_indices[i]](0),
                                         m_particles[particle_indices[i]](1),
                                         m_particles[particle_indices[i]](2));
        } else {
          dofs.segment<4>(i * 4) = m_particles[particle_indices[i]];
        }
      }

      ElasticStrand* strand =
          new ElasticStrand(dofs, area_dofs, m_strand_parameters[paramsIndex],
                            m_collision_parameters[paramsIndex], controller);

      strand->setGlobalIndex(numstrands);
      m_strands.push_back(strand);

      RodData* rd = new RodData(*strand, *controller, num_components);
      m_rodDatum.push_back(rd);

      m_local_to_global.push_back(particle_indices);

      ++numstrands;
      ++linenumber;
    }

    ifs.close();
    ++i_obj_strand;
  }
}

void dump_binary_checkpoint_subprog(DumpDataBinary* data) {
  std::ofstream os_data(data->data_fn.c_str(), std::ios::binary);

  // write general info
  os_data.write((char*)&(data->m_time), sizeof(Scalar));
  os_data.write((char*)&(data->current_frame), sizeof(int));
  os_data.write((char*)&(data->current_checkpoint), sizeof(int));
  os_data.write((char*)&(data->num_components), sizeof(int));

  // write fluid
  const int num_particles = data->m_x.size() / 3;
  os_data.write((const char*)&num_particles, sizeof(int));

  if (num_particles) {
    os_data.write((char*)data->m_x.data(), num_particles * 3 * sizeof(Scalar));
    os_data.write((char*)data->m_v.data(), num_particles * 3 * sizeof(Scalar));
    os_data.write((char*)data->m_m.data(), num_particles * 3 * sizeof(Scalar));
    os_data.write((char*)data->m_radius.data(), num_particles * sizeof(Scalar));
    os_data.write((char*)data->m_J.data(), num_particles * sizeof(Scalar));
    os_data.write((char*)data->m_vol.data(), num_particles * sizeof(Scalar));
    os_data.write((char*)data->m_rest_vol.data(),
                  num_particles * sizeof(Scalar));
    os_data.write((char*)&(data->m_particle_group[0]),
                  num_particles * sizeof(int));
    os_data.write((char*)&(data->m_classifier[0]),
                  num_particles * sizeof(ParticleClassifier));
    os_data.write((char*)data->m_weakened.data(),
                  num_particles * sizeof(unsigned char));
    os_data.write((char*)data->m_components.data(),
                  num_particles * data->num_components * sizeof(Scalar));
    os_data.write((char*)data->m_proj_func.data(),
                  num_particles * 3 * sizeof(Scalar));

    os_data.write((char*)data->m_Fe.data(), num_particles * 9 * sizeof(Scalar));
    os_data.write((char*)data->m_b.data(), num_particles * 9 * sizeof(Scalar));
    os_data.write((char*)data->m_B.data(), num_particles * 9 * sizeof(Scalar));
    os_data.write((char*)data->m_b_trial.data(),
                  num_particles * 9 * sizeof(Scalar));
  }

  // write hairs
  const int num_strands = data->m_currentDOFs.size();
  os_data.write((char*)&num_strands, sizeof(int));

  for (int i = 0; i < num_strands; ++i) {
    const int num_verts = data->m_currentAreaDOFs[i].size();
    os_data.write((char*)&num_verts, sizeof(int));

    os_data.write((char*)(data->m_currentDOFs[i].data()),
                  data->m_currentDOFs[i].size() * sizeof(Scalar));
    os_data.write((char*)(data->m_currentAreaDOFs[i].data()),
                  data->m_currentAreaDOFs[i].size() * sizeof(Scalar));
    os_data.write((char*)(data->m_velocities[i].data()),
                  data->m_velocities[i].size() * sizeof(Scalar));
    os_data.write((char*)(data->m_flow_velocities[i].data()),
                  data->m_flow_velocities[i].size() * sizeof(Scalar));
    os_data.write((char*)(data->m_flow_strain[i].data()),
                  data->m_flow_strain[i].size() * sizeof(Scalar));
    os_data.write((char*)(data->m_flow_components[i].data()),
                  data->m_flow_components[i].size() * sizeof(Scalar));
    os_data.write((char*)(data->m_flow_reservoirs[i].data()),
                  data->m_flow_reservoirs[i].size() * sizeof(Scalar));
  }

  const int num_goals = data->m_strand_goals.size();
  os_data.write((char*)&num_goals, sizeof(int));
  if (num_goals) {
    os_data.write((char*)&(data->m_strand_goals[0]),
                  sizeof(std::pair<Vec2i, Vec4x>) * num_goals);
  }

  // write groups
  const int num_groups = data->m_group_pos.size();
  os_data.write((char*)&num_groups, sizeof(int));

  if (num_groups) {
    os_data.write((char*)&(data->m_group_pos[0]), sizeof(Vec3x) * num_groups);
    os_data.write((char*)&(data->m_group_scale[0]), sizeof(Vec3x) * num_groups);
    os_data.write((char*)&(data->m_group_rot[0]),
                  sizeof(Eigen::Quaternion<Scalar>) * num_groups);
  }

  // write fields
  const int num_fields = data->m_field_i_frames.size();
  os_data.write((char*)&num_fields, sizeof(int));

  if (num_fields) {
    os_data.write((char*)&(data->m_field_i_frames[0]),
                  sizeof(int) * num_fields);
    os_data.write((char*)&(data->m_field_center[0]),
                  sizeof(Vec3x) * num_fields);
    os_data.write((char*)&(data->m_field_rot[0]),
                  sizeof(Eigen::Quaternion<Scalar>) * num_fields);
    os_data.write((char*)&(data->m_field_future_scale[0]),
                  sizeof(Vec3x) * num_fields);
  }

  for (int i = 0; i < num_fields; ++i) {
    const int num_params = data->m_field_parameter[i].size();
    os_data.write((char*)&(num_params), sizeof(int));
    os_data.write((char*)(data->m_field_parameter[i].data()),
                  data->m_field_parameter[i].size() * sizeof(Scalar));
  }

  // write mesh
  const int num_meshes = data->m_mesh_current_frames.size();
  os_data.write((char*)&num_meshes, sizeof(int));

  for (int i = 0; i < num_meshes; ++i) {
    const int num_verts = data->m_previous_mesh_vertices[i].size();
    os_data.write((char*)&num_verts, sizeof(int));

    os_data.write((char*)&(data->m_mesh_current_frames[i]), sizeof(int));
    os_data.write((char*)&(data->m_previous_mesh_vertices[i][0]),
                  sizeof(Vec3x) * data->m_previous_mesh_vertices[i].size());
    os_data.write((char*)&(data->m_current_mesh_vertices[i][0]),
                  sizeof(Vec3x) * data->m_current_mesh_vertices[i].size());
    os_data.write((char*)&(data->m_next_mesh_vertices[i][0]),
                  sizeof(Vec3x) * data->m_next_mesh_vertices[i].size());
    os_data.write((char*)&(data->m_current_mesh_displacement[i][0]),
                  sizeof(Vec3x) * data->m_current_mesh_displacement[i].size());
  }

  os_data.flush();
  os_data.close();
  delete data;
}

void XMLReader::dumpBinaryCheckpoint(std::string outputdirectory,
                                     int current_frame, int current_checkpoint,
                                     int file_width) const {
  std::stringstream name;
  name << std::setfill('0');
  name << outputdirectory << "/checkpoint_" << std::setw(file_width)
       << current_checkpoint << ".bin";

  DumpDataBinary* data = new DumpDataBinary;

  data->data_fn = name.str().c_str();
  data->m_time = m_stepper->getTime();
  data->current_frame = current_frame;
  data->current_checkpoint = current_checkpoint;
  data->num_components = m_num_components;

  // fill fluid variables
  if (m_fluidScripting_controllers.size() && m_fluidScripting_controllers[0]) {
    std::shared_ptr<FluidScriptingController> controller =
        m_fluidScripting_controllers[0];

    data->m_x = controller->getX();
    data->m_v = controller->getV();
    data->m_m = controller->getM();
    data->m_radius = controller->getRadius();
    data->m_J = controller->getJ();
    data->m_vol = controller->getVol();
    data->m_rest_vol = controller->getRestVol();
    data->m_particle_group = controller->getParticleGroup();
    data->m_classifier = controller->getClassifier();
    data->m_weakened = controller->getWeakened();
    data->m_components = controller->getComponents();
    data->m_proj_func = controller->getProjFunc();

    data->m_Fe = controller->getFe();
    data->m_b = controller->getb();
    data->m_B = controller->getB();
    data->m_b_trial = controller->getbtrial();
  }

  // fill rod variables
  const int num_rods = m_rodDatum.size();
  data->m_currentDOFs.resize(num_rods);
  data->m_currentAreaDOFs.resize(num_rods);
  data->m_velocities.resize(num_rods);
  data->m_flow_velocities.resize(num_rods);
  data->m_flow_strain.resize(num_rods);
  data->m_flow_components.resize(num_rods);
  data->m_flow_reservoirs.resize(num_rods);

  for (int i = 0; i < num_rods; ++i) {
    ElasticStrand& strand = m_rodDatum[i]->getStrand();
    data->m_currentDOFs[i] = strand.getCurrentDegreesOfFreedom();
    data->m_currentAreaDOFs[i] = strand.getCurrentAreaDegreesOfFreedom();
    data->m_velocities[i] = strand.getStepper()->velocities();
    data->m_flow_velocities[i] = strand.getStepper()->getCurrentFlowVelocity();
    data->m_flow_strain[i] = strand.getStepper()->flowStrain();
    data->m_flow_components[i] = strand.getStepper()->flowComponents();
    data->m_flow_reservoirs[i] = strand.getReservoir();

    DOFScriptingController* controller =
        strand.dynamics().getScriptingController();

    const int num_vtx = strand.getNumVertices();
    for (int j = 0; j < num_vtx; ++j) {
      if (controller->isVertexGoaled(j) || controller->isThetaGoaled(j)) {
        Vec3x v = controller->getVertexGoal(j);
        Scalar t = controller->getThetaGoal(j);
        data->m_strand_goals.push_back(
            std::pair<Vec2i, Vec4x>(Vec2i(i, j), Vec4x(v(0), v(1), v(2), t)));
      }
    }
  }

  // fill group variables;
  data->m_group_pos = m_group_pos;
  data->m_group_scale = m_group_scale;
  data->m_group_rot = m_group_rot;

  const int num_dfos =
      m_fields.size() + m_sources.size() + m_terminators.size();
  data->m_field_i_frames.reserve(num_dfos);
  data->m_field_center.reserve(num_dfos);
  data->m_field_parameter.reserve(num_dfos);
  data->m_field_rot.reserve(num_dfos);
  data->m_field_future_scale.reserve(num_dfos);

  // fill field variables
  for (const DistanceFieldObject& dfo : m_fields) {
    data->m_field_i_frames.push_back(dfo.i_frame);
    data->m_field_center.push_back(dfo.center);
    data->m_field_parameter.push_back(dfo.parameter);
    data->m_field_rot.push_back(dfo.rot);
    data->m_field_future_scale.push_back(dfo.future_scale);
  }

  for (const DistanceFieldObject& dfo : m_sources) {
    data->m_field_i_frames.push_back(dfo.i_frame);
    data->m_field_center.push_back(dfo.center);
    data->m_field_parameter.push_back(dfo.parameter);
    data->m_field_rot.push_back(dfo.rot);
    data->m_field_future_scale.push_back(dfo.future_scale);
  }

  for (const DistanceFieldObject& dfo : m_terminators) {
    data->m_field_i_frames.push_back(dfo.i_frame);
    data->m_field_center.push_back(dfo.center);
    data->m_field_parameter.push_back(dfo.parameter);
    data->m_field_rot.push_back(dfo.rot);
    data->m_field_future_scale.push_back(dfo.future_scale);
  }

  const int num_meshes = m_meshScripting_controllers.size();
  // fill mesh variables
  data->m_mesh_current_frames.resize(num_meshes);
  data->m_previous_mesh_vertices.resize(num_meshes);
  data->m_current_mesh_vertices.resize(num_meshes);
  data->m_next_mesh_vertices.resize(num_meshes);
  data->m_current_mesh_displacement.resize(num_meshes);

  for (int i = 0; i < num_meshes; ++i) {
    data->m_mesh_current_frames[i] =
        m_meshScripting_controllers[i]->getIFrame();
    data->m_previous_mesh_vertices[i] =
        m_meshScripting_controllers[i]->getPreviousMesh()->m_vertices;
    data->m_current_mesh_vertices[i] =
        m_meshScripting_controllers[i]->getCurrentMesh()->m_vertices;
    data->m_next_mesh_vertices[i] =
        m_meshScripting_controllers[i]->getNextMesh()->m_vertices;
    data->m_current_mesh_displacement[i] =
        m_meshScripting_controllers[i]->getCurrentMesh()->m_displacements;
  }

  std::thread(std::bind(dump_binary_checkpoint_subprog, data)).detach();
}

void XMLReader::apply_global_scaling(int group_idx, const Vec3x& mul) {
  for (int i : m_group_fields[group_idx]) {
    m_fields[i].apply_global_scaling(mul);
  }

  for (int i : m_group_sources[group_idx]) {
    m_sources[i].apply_global_scaling(mul);
  }
}

void XMLReader::apply_local_scaling(int group_idx, const Vec3x& mul) {
  for (int i : m_group_fields[group_idx]) {
    m_fields[i].apply_local_scaling(mul);
  }

  for (int i : m_group_sources[group_idx]) {
    m_sources[i].apply_local_scaling(mul);
  }
}

void XMLReader::apply_global_rotation(int group_idx,
                                      const Eigen::Quaternion<double>& rot) {
  for (int i : m_group_fields[group_idx]) {
    m_fields[i].apply_global_rotation(rot);
  }

  for (int i : m_group_sources[group_idx]) {
    m_sources[i].apply_global_rotation(rot);
  }
}

void XMLReader::apply_local_rotation(int group_idx,
                                     const Eigen::Quaternion<double>& rot) {
  for (int i : m_group_fields[group_idx]) {
    m_fields[i].apply_local_rotation(rot);
  }

  for (int i : m_group_sources[group_idx]) {
    m_sources[i].apply_global_rotation(rot);
  }
}

void XMLReader::apply_translation(int group_idx, const Vec3x& t) {
  for (int i : m_group_fields[group_idx]) {
    m_fields[i].apply_translation(t);
  }

  for (int i : m_group_sources[group_idx]) {
    m_sources[i].apply_translation(t);
  }
}

void XMLReader::apply_switch(int group_idx, unsigned mask) {
  const std::vector<std::pair<int, int> >& group = m_groups[group_idx];
  for (const std::pair<int, int>& p : group) {
    RodData* rd = m_rodDatum[p.first];

    if (mask & 1U) {
      if (rd->getDofController().isVertexGoaled(p.second) ||
          rd->getDofController().isVertexFreezed(p.second)) {
        if (m_simulation_params.m_useSoftAttachConstraints) {
          rd->getDofController().ungoalVertices(p.second);
        } else {
          rd->getDofController().unfreezeVertices(p.second);
        }

        m_fixed[m_local_to_global[p.first][p.second]] &= 0xFEU;
      } else {
        if (m_simulation_params.m_useSoftAttachConstraints) {
          rd->getDofController().setVertexGoal(
              p.second, rd->getStrand().getVertex(p.second));
        } else {
          rd->getDofController().freezeVertices(p.second);
        }

        m_fixed[m_local_to_global[p.first][p.second]] |= 0x1U;
      }
    }

    if (mask & 2U) {
      if (rd->getDofController().isThetaGoaled(p.second) ||
          rd->getDofController().isThetaFreezed(p.second)) {
        if (m_simulation_params.m_useSoftAttachConstraints) {
          rd->getDofController().ungoalTheta(p.second);
        } else {
          rd->getDofController().unfreezeTheta(p.second);
        }

        m_fixed[m_local_to_global[p.first][p.second]] &= 0xFDU;
      } else {
        if (m_simulation_params.m_useSoftAttachConstraints) {
          rd->getDofController().setThetaGoal(
              p.second, rd->getStrand().getTheta(p.second));
        } else {
          rd->getDofController().freezeTheta(p.second);
        }

        m_fixed[m_local_to_global[p.first][p.second]] |= 0x2U;
      }
    }
  }
}

void XMLReader::projectConstraint(
    const std::vector<ImplicitStepper*>& steppers) {}

double Script::getNextVelocity(const double& dt, const double& current_time) {
  if (func == Script::COSINE) {
    return cosine_ease_function(current_time + dt, start, end,
                                start + ease_start, end - ease_end, amplitude,
                                frequency);
  } else if (func == Script::CUBIC) {
    return cubic_ease_function(current_time + dt, start, end,
                               start + ease_start, end - ease_end, v(3));
  } else if (func == Script::WENO) {
    const double vel = weno_ease_function(current_time + dt, dt, start, end,
                                          base_dt, base_pos, base_vertices);
    base_pos += vel * dt;
    return vel;
  }

  return 0.0;
}

void Script::stepScript(const double& dt, const double& current_time) {
  if (!m_scene) return;

  if (current_time >= start && current_time + dt <= end) {
    switch (type) {
      case Script::SCALE: {
        if (func != Script::CUBIC) {
          std::cout << "[error: non-cubic scaling is unsupported!]"
                    << std::endl;
          break;
        }

        Vec3x& scale = m_scene->m_group_scale[group_index];
        m_scene->m_group_prev_scale[group_index] = scale;

        Vec3x vel = Vec3x(
            cubic_ease_function(current_time + dt, start, end,
                                start + ease_start, end - ease_end, log(v(0))),
            cubic_ease_function(current_time + dt, start, end,
                                start + ease_start, end - ease_end, log(v(1))),
            cubic_ease_function(current_time + dt, start, end,
                                start + ease_start, end - ease_end, log(v(2))));

        Vec3x multipliers =
            Vec3x(exp(vel(0) * dt), exp(vel(1) * dt), exp(vel(2) * dt));

        scale.array() *= multipliers.array();
        if (transform_global) {
          m_scene->apply_global_scaling(group_index, multipliers);

          Vec3x& trans = m_scene->m_group_pos[group_index];
          m_scene->m_group_prev_pos[group_index] = trans;

          trans.array() *= multipliers.array();
        } else {
          m_scene->apply_local_scaling(group_index, multipliers);
        }

        break;
      }
      case Script::SWITCH: {
        if (!used) {
          m_scene->apply_switch(group_index, switch_mask);
          used = true;
        }
        break;
      }
      case Script::TRANSLATE: {
        Vec3x vec(v(0), v(1), v(2));
        double dx = vec.norm();
        v(3) = dx;

        Vec3x& trans = m_scene->m_group_pos[group_index];
        m_scene->m_group_prev_pos[group_index] = trans;
        if (dx > 0.0) {
          VecXx nv = vec / dx;

          double vel = getNextVelocity(dt, current_time);
          trans += nv * vel * dt;
          m_scene->apply_translation(group_index, nv * vel * dt);
        }

        break;
      }
      case Script::ROTATE: {
        double vel = getNextVelocity(dt, current_time);
        double da = vel * dt;
        Vec3x axis(v(0), v(1), v(2));

        Eigen::AngleAxis<double> rot(da, axis);
        Eigen::Quaternion<double> qrot(rot);

        Eigen::Quaternion<double>& prot = m_scene->m_group_rot[group_index];
        m_scene->m_group_prev_rot[group_index] = prot;

        if (transform_with_origin) {
          Vec3x& trans = m_scene->m_group_pos[group_index];
          m_scene->m_group_prev_pos[group_index] = trans;
        }

        if (transform_global) {
          prot = qrot * prot;
          m_scene->apply_global_rotation(group_index, qrot);
        } else {
          prot *= qrot;
          m_scene->apply_local_rotation(group_index, qrot);
        }

        if (transform_with_origin) {
          const Eigen::Quaternion<double> q_diff =
              prot * m_scene->m_group_prev_rot[group_index].inverse();
          Vec3x& trans = m_scene->m_group_pos[group_index];
          Vec3x dir_vec = trans - origin;

          Vec3x rot_dir_vec = q_diff * dir_vec;
          Vec3x vec = rot_dir_vec - dir_vec;
          trans += vec;
          m_scene->apply_translation(group_index, vec);
        }
        break;
      }
      case Script::SWIRL: {
        double vel = getNextVelocity(dt, current_time);
        double da = vel * dt;
        Vec3x axis(v(0), v(1), v(2));

        Eigen::AngleAxis<double> rot(da, axis);
        Eigen::Quaternion<double> qrot(rot);

        Eigen::Quaternion<double>& prot = m_scene->m_group_rot[group_index];
        m_scene->m_group_prev_rot[group_index] = prot;

        prot = qrot * prot;

        Vec3x& trans = m_scene->m_group_pos[group_index];
        m_scene->m_group_prev_pos[group_index] = trans;

        const Eigen::Quaternion<double> q_diff =
            prot * m_scene->m_group_prev_rot[group_index].inverse();
        Vec3x dir_vec = trans - origin;

        Vec3x rot_dir_vec = q_diff * dir_vec;
        Vec3x vec = rot_dir_vec - dir_vec;
        trans += vec;
        m_scene->apply_translation(group_index, vec);
        break;
      }
      default:
        std::cout << "UNIMPLEMENTED SCRIPT TYPE [" << type << "]!" << std::endl;
        break;
    }
  }
}
