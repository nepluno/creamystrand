cmake_minimum_required( VERSION 2.8 )

project( StrandSim C CXX )


###############################################################################
### Configure and check the build setting options

# If the user provided a CMAKE_BUILD_TYPE, create a cache variable
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

if (APPLE)
	option(USE_XCODE_PATCH "This is a option for XCode patch on MAC" ON)
endif()

###############################################################################
### Allow the user to set a preferred BLAS and LAPACK implementation

set( SUGGESTED_BLAS_LAPACK "" CACHE STRING "The desired LAPACK implementation (e.g. Generic, Apple, ACML, Intel10_64lp, Intel10_64lp_seq, ..." )


###############################################################################
### Set flags used across all builds (TODO: Check compiler types here!)

# Append an extra space if flags are already set
if( CMAKE_CXX_FLAGS )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} " )
endif( CMAKE_CXX_FLAGS )

# Enable c++0x support
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DEIGEN_DONT_ALIGN -DEIGEN_DONT_ALIGN_STATICALLY" )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC" )

###############################################################################
### Set preprocessor defines used in all builds

if (WIN32)
  add_definitions("-DBOOST_ALL_NO_LIB")
endif()

###############################################################################
### Set preprocessor defines used in debug builds

if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  add_definitions( -DEIGEN_DONT_ALIGN ) # Also does -DEIGEN_DONT_VECTORIZE
endif( CMAKE_BUILD_TYPE STREQUAL "Debug" )


###############################################################################
### Set preprocessor defines used in release builds

#if( CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" OR CMAKE_BUILD_TYPE STREQUAL "MinSizeRel" )
#  add_definitions(-Dother_flags_go_here)
#endif( CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" OR CMAKE_BUILD_TYPE STREQUAL "MinSizeRel" )


###############################################################################
### Set paths, etc

# Create a default location for custom cmake modules
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake )

# Allow projects to reference source relative to the root of the project
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

###############################################################################
### Add libraries, plugins

add_subdirectory( bogus )
add_subdirectory( StrandSim )
add_subdirectory( Apps/StrandSimulator )
