/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#ifndef __RENDERING_UTILITIES_H__
#define __RENDERING_UTILITIES_H__

#ifdef WIN32
#include <Windows.h>
#endif

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <list>
#include <iostream>
#include <cstdio>

#include "../../StrandSim/Core/Definitions.hh"
#include "../../StrandSim/Utils/StringUtilities.hh"
#include "../../StrandSim/Utils/MathUtilities.hh"
#include "JET.hh"

namespace renderingutils
{
    using namespace strandsim;
	// False => error
	bool checkGLErrors();
	
	inline Vec3x interpolateColor(const Scalar& x, const Scalar xmin = 0.0, const Scalar xmax = 1.0)
	{
		Scalar dm = (xmax - xmin);
		
		Scalar a;
		if(dm == 0.0) a = x;
		else a = (x - xmin) / dm * (Scalar)(jetmapping_size - 1);
		
		int isel = std::max(std::min((int) a, jetmapping_size - 1), 0);
		int inext = (isel + 1) % (jetmapping_size);
		Scalar fraca = std::max(std::min(a - (Scalar) isel, 1.0), 0.0);
		
		return lerp(jetmapping_real[isel], jetmapping_real[inext], fraca);
	}
}

#endif
