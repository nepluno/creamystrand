/**
 * \copyright 2010 Breannan Smith
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#include "RenderingUtilities.hh"

namespace renderingutils
{
	bool checkGLErrors()
	{
//#ifdef RENDER_ENABLED
//        GLenum errCode;
//        const GLubyte *errString;
//
//        if ((errCode = glGetError()) != GL_NO_ERROR)
//        {
//            errString = gluErrorString(errCode);
//            std::cout << outputmod::startred << "OpenGL Error:" << outputmod::endred << std::flush;
//            fprintf(stderr, " %s\n", errString);
//            return false;
//        }
//#endif
		return true;
	}
}
