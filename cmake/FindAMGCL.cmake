# - Try to find AMGCL lib
#
# Once done this will define
#
#  AMGCL_FOUND - system has eigen lib with correct version
#  AMGCL_INCLUDE_DIR - the eigen include directory

find_path(AMGCL_INCLUDE_DIR NAMES amgcl/amg.hpp
    PATHS
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_INSTALL_PREFIX}/include
    ${CMAKE_SOURCE_DIR}/External
    ${CMAKE_INSTALL_PREFIX}/External
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AMGCL DEFAULT_MSG AMGCL_INCLUDE_DIR)
mark_as_advanced(AMGCL_INCLUDE_DIR)
