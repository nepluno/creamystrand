# - Find Visualization Library
# Find the Visualization Library headers and libraries
#
# VL_INCLUDE_DIR - include path for headers
# VL_LIBRARIES   - libraries to include when linking
# VL_FOUND       - True if VisualizationLibrary is found

if (VL_INCLUDE_DIR AND VL_LIBRARIES)
  # already in cache, be silent
  set (VL_FIND_QUIETLY TRUE)
endif (VL_INCLUDE_DIR AND VL_LIBRARIES)

# find the headers
find_path (VL_INCLUDE_PATH
  vl/VisualizationLibrary.hpp
  HINTS ENV VisualizationLibrary_INC_DIR
  DOC "Visualization Library include directory (http://www.visualizationlibrary.com/)"
  )

find_library (VL_VisualizationLibrary_LIBRARY
  NAMES VisualizationLibrary
  HINTS ENV VisualizationLibrary_LIB_DIR
  )

find_library (VL_VLGLUT_LIBRARY
  NAMES VLGLUT
  HINTS ENV VisualizationLibrary_LIB_DIR
  )

find_library (VL_FreeType_LIBRARY
  NAMES FreeType
  HINTS ENV VisualizationLibrary_LIB_DIR
  )

find_library (VL_JPG_LIBRARY
  NAMES JPG
  HINTS ENV VisualizationLibrary_LIB_DIR
  )

find_library (VL_PNG_LIBRARY
  NAMES PNG
  HINTS ENV VisualizationLibrary_LIB_DIR
  )

find_library (VL_TIFF_LIBRARY
  NAMES TIFF
  HINTS ENV VisualizationLibrary_LIB_DIR
  )

find_library (VL_ZLib_LIBRARY
  NAMES ZLib
  HINTS ENV VisualizationLibrary_LIB_DIR
  )

# handle the QUIETLY and REQUIRED arguments and set
# VL_FOUND to TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (VL "VisualizationLibrary (http://www.visualizationlibrary.com/) not found. Set VisualizationLibarary_INC_DIR and VisualizationLibrary_LIB_DIR." VL_INCLUDE_PATH VL_VisualizationLibrary_LIBRARY VL_VLGLUT_LIBRARY VL_FreeType_LIBRARY VL_PNG_LIBRARY VL_JPG_LIBRARY VL_TIFF_LIBRARY VL_ZLib_LIBRARY)

if (VL_FOUND)
  set (VL_INCLUDE_DIR ${VL_INCLUDE_PATH})
  set (VL_LIBRARIES ${VL_VLGLUT_LIBRARY} ${VL_FreeType_LIBRARY} ${VL_JPG_LIBRARY} ${VL_PNG_LIBRARY} ${VL_TIFF_LIBRARY} ${VL_ZLib_LIBRARY} ${VL_VisualizationLibrary_LIBRARY})
endif (VL_FOUND)
