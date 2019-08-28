/**
 * \copyright 2010 Breannan Smith
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TRIANGLEMESHRENDERER_HH
#define TRIANGLEMESHRENDERER_HH


#include "../../../StrandSim/Collision/TriangularMesh.hh"
#include "../../StrandSim/Render/RenderBase.hh"

namespace strandsim {
    
    class FluidScriptingController;
    /** Class that implements OpenGL rendering for triangle meshes. */
    class TriangleMeshRenderer : public RenderBase
    {
    public:
        
        enum DrawMode { DBG, FLAT, NONE };
        
        explicit TriangleMeshRenderer( TriangularMesh& mesh, int num_components );
        
        void render();
        
        DrawMode getMode() const { return m_mode; }
        void setMode(DrawMode mode) { m_mode = mode; }
        
        virtual Vec3d calculateObjectCenter();
        virtual Scalar calculateObjectBoundingRadius(const Vec3d& center);
        
        void setFluidController( const std::shared_ptr< FluidScriptingController >& controller );
        
    protected:
        TriangularMesh& m_mesh;
        DrawMode m_mode;
        
        std::vector< Vec3x > m_component_colors;
        
        std::shared_ptr< FluidScriptingController > m_controller;
    };
    
}

#endif 
