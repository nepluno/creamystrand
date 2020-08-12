/**
 * \copyright 2010 Breannan Smith
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MESHSCRIPTINGCONTROLLER_HH
#define MESHSCRIPTINGCONTROLLER_HH

#include "../Forces/LevelSetFwd.hh"
#include "ScriptingController.hh"

namespace strandsim {

class TriangularMesh;
class TriangularMeshFlow;

class MeshScriptingController
    : public ScriptingController,
      public std::enable_shared_from_this<MeshScriptingController> {
 public:
  MeshScriptingController(Scalar time, Scalar dt);
  virtual ~MeshScriptingController();

  virtual void createInitialLevelSet() = 0;

  virtual bool execute(bool updateLevelSet) = 0;
  virtual bool execute();

  virtual strandsim::LevelSet* currentLevelSet() = 0;
  virtual bool hasLevelSet() const = 0;

  virtual std::shared_ptr<TriangularMesh> getCurrentMesh() = 0;
  virtual std::shared_ptr<TriangularMesh> getPreviousMesh() = 0;
  virtual std::shared_ptr<TriangularMesh> getNextMesh() = 0;
  virtual int getIFrame() = 0;
  virtual void setIFrame(int iframe) = 0;
  // virtual const TriangularMesh* getCurrentMesh() const = 0;
  virtual const std::vector<bool>& getEnabledVertices() const = 0;

  virtual double getDefaultFrictionCoefficient() const = 0;
  virtual const std::vector<double>& getFrictionCoefficients() const = 0;
  virtual double getLevelSetForceThickness() const = 0;
  virtual double getLevelSetForceStrength() const = 0;
  virtual std::shared_ptr<TriangularMeshFlow> getMeshFlow() const = 0;

  virtual bool updateMeshNormalArea() = 0;

  // Virtual functions that should be inherited for meshes that are collidable
  // on both sides
  virtual bool collideOnBothSides() const { return false; }

  // Whether the controller is sure of the side of which a vertex should stay.
  // Returns: param atPreviousStep ; if true, do not tka einto account guesses
  // from the current step
  //  1 : normal is known and correspond to counter-clockwise vector product of
  //  vertices 0 : normal is unknown
  // -1 : normal is known and correspond to clockwise vector product of vertices
  virtual short knowsNormalSign(bool atPreviousStep, unsigned faceIndex,
                                unsigned rodIndex, unsigned vertex) {
    return 1;
  }
  virtual void setNormalSign(short sign, float offset, unsigned faceIndex,
                             unsigned rodIndex, unsigned vertex) {}
};

}  // namespace strandsim

#endif
