/**
 * \copyright 2010 Breannan Smith
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#ifndef FLUIDSRENDERER_HH
#define FLUIDSRENDERER_HH


#include "../../StrandSim/Render/RenderBase.hh"
#include "../../StrandSim/Dynamic/FluidScriptingController.hh"
#include "../../StrandSim/Utils/ThreadUtils.hh"

namespace strandsim {
    
    /** Class that implements OpenGL rendering for triangle meshes. */
    class FluidsRenderer : public RenderBase
    {
    public:
        
        enum DrawMode { DBG, FANCY };
        
        explicit FluidsRenderer( const std::shared_ptr<FluidScriptingController>& pcontroller, MutexType& particle_mutex, MutexType& grid_mutex, FluidsRenderer::DrawMode dm = DBG );
        
        void render();
        
        DrawMode getMode() const { return m_mode; }
        void setMode(DrawMode mode) { m_mode = mode; }
        
    protected:
        const std::shared_ptr<FluidScriptingController> m_pcontroller;
        DrawMode m_mode;
        MutexType& m_particle_mutex;
        MutexType& m_grid_mutex;
        
        std::vector< Vec3x > m_component_colors;
    };
    
}

#endif 
