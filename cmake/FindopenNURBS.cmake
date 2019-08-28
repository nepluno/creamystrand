# - Find openNURBS
# Find the openNURBS headers and libraries
#
# OPENNURBS_INCLUDE_DIR - include path for headers
# OPENNURBS_LIBRARIES   - libraries to include when linking
# OPENNURBS_FOUND       - True if openNURBS is found

if( OPENNURBS_INCLUDE_DIR AND OPENNURBS_LIBRARIES )
  # already in cache, be silent
  set( OPENNURBS_FIND_QUIETLY TRUE )
endif( OPENNURBS_INCLUDE_DIR AND OPENNURBS_LIBRARIES )

# find the headers
find_path( OPENNURBS_INCLUDE_PATH
  openNURBS/opennurbs.h
  HINTS ${CMAKE_SOURCE_DIR}/External
  )

find_library( OPENNURBS_LIBRARY
  NAMES openNURBS
  HINTS ${CMAKE_SOURCE_DIR}/External/openNURBS/lib
  )

# handle the QUIETLY and REQUIRED arguments and set
include ( FindPackageHandleStandardArgs )
find_package_handle_standard_args( openNURBS "openNURBS (http://www.opennurbs.org/) not found." OPENNURBS_INCLUDE_PATH OPENNURBS_LIBRARY )

if( OPENNURBS_FOUND )
  set( OPENNURBS_INCLUDE_DIR ${OPENNURBS_INCLUDE_PATH} )
  set( OPENNURBS_LIBRARIES ${OPENNURBS_LIBRARY} )
endif( OPENNURBS_FOUND )

mark_as_advanced( OPENNURBS_INCLUDE_PATH )
mark_as_advanced( OPENNURBS_LIBRARY )
