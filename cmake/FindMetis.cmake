# - Find metis
# Find the metis headers and libraries
#
# Metis_INCLUDE_DIR - include path for headers
# Metis_LIBRARIES   - libraries to include when linking
# Metis_FOUND       - True if metis is found

if( Metis_INCLUDE_DIR AND Metis_LIBRARIES )
  # already in cache, be silent
  set( Metis_FIND_QUIETLY TRUE )
endif( Metis_INCLUDE_DIR AND Metis_LIBRARIES )

# find the headers
find_path( Metis_INCLUDE_PATH
  metis.h
  HINTS ${CMAKE_SOURCE_DIR}/External/metis $ENV{METIS_PATH}/include
  PATH_SUFFIXES metis
  )

find_library( Metis_LIBRARY
  NAMES metis
  HINTS ${CMAKE_SOURCE_DIR}/External/metis $ENV{METIS_PATH}/lib
  )

# handle the QUIETLY and REQUIRED arguments and set
include ( FindPackageHandleStandardArgs )
find_package_handle_standard_args( Metis "Metis (http://www.cs.umn.edu/~metis/) not found." Metis_INCLUDE_PATH Metis_LIBRARY )

if( METIS_FOUND )
  set( Metis_INCLUDE_DIR ${Metis_INCLUDE_PATH} )
  set( Metis_LIBRARIES ${Metis_LIBRARY} )
endif( METIS_FOUND )

mark_as_advanced( Metis_INCLUDE_PATH )
mark_as_advanced( Metis_LIBRARY )
