# - Find Approximate Nearest Neighbors
# Find the Approximate Nearest Neighbors headers and libraries
#
# If your ANN installation is not in a standard location (e.g. /usr/include)
# set the environment variable ANN_PATH to point to your install.
#
# ANN_INCLUDE_DIR - include path for headers
# ANN_LIBRARIES   - libraries to include when linking
# ANN_FOUND       - True if ANN is found

if( ANN_INCLUDE_DIR AND ANN_LIBRARIES )
  # already in cache, be silent
  set( ANN_FIND_QUIETLY TRUE )
endif( ANN_INCLUDE_DIR AND ANN_LIBRARIES )

# find the headers
find_path( ANN_INCLUDE_PATH
  ANN/ANN.h
  HINTS ${CMAKE_SOURCE_DIR}/External/ann/include $ENV{ANN_PATH}/include
  )

find_library( ANN_LIBRARY
  NAMES ANN
  HINTS ${CMAKE_SOURCE_DIR}/External/ann/lib $ENV{ANN_PATH}/lib
  )

if(NOT ANN_LIBRARY )
find_library( ANN_LIBRARY
  NAMES ann
  HINTS ${CMAKE_SOURCE_DIR}/External/ann/lib $ENV{ANN_PATH}/lib
  )
endif()

# handle the QUIETLY and REQUIRED arguments and set
include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( ANN "ANN (http://www.cs.umd.edu/~mount/ANN/) not found." ANN_INCLUDE_PATH ANN_LIBRARY )
#set( ANN_INCLUDE_DIR /vol/bob/built/linux64/L43/ext/ANN/tags/ANN-1.0/include )
#set( ANN_LIBRARIES /vol/bob/built/linux64/L43/ext/ANN/tags/ANN-1.0/lib )

if( ANN_FOUND )
  set( ANN_INCLUDE_DIR ${ANN_INCLUDE_PATH} )
  set( ANN_LIBRARIES ${ANN_LIBRARY} )
endif( ANN_FOUND )

# Clear internal variables
mark_as_advanced( ANN_INCLUDE_PATH )
mark_as_advanced( ANN_LIBRARY )
