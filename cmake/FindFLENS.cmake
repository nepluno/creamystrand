# - Find FLENS
# Find the FLENS includes and library
#
# FLENS_INCLUDE_DIR - where to find the FLENS headers
# FLENS_LIBRARIES   - list of libraries when using FLENS
# FLENS_FOUND       - True if FLENS is found

if (FLENS_INCLUDE_DIR AND FLENS_LIBRARIES)
  # already in cache, be silent
  set (FLENS_FIND_QUIETLY TRUE)
endif (FLENS_INCLUDE_DIR AND FLENS_LIBRARIES)

# find the headers
find_path (FLENS_INCLUDE_PATH
  flens/flens.h
  HINTS FLENS_INC_DIR ENV FLENS_INC_DIR
  )

# find the libraries
set (FLENS_NAMES flens)
find_library (FLENS_LIBRARY NAMES
  NAMES ${FLENS_NAMES}
  HINTS FLENS_LIB_DIR ENV FLENS_LIB_DIR
  )

# handle the QUIETLY and REQUIRED arguments and set FLENS_FOUND to
# TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FLENS "FLENS (http://flens.sourceforge.net/) could not be found. Set FLENS_INC_DIR and FLENS_LIB_DIR." FLENS_INCLUDE_PATH FLENS_LIBRARY)

if (FLENS_FOUND)
  set (FLENS_INCLUDE_DIR ${FLENS_INCLUDE_PATH})
  set (FLENS_LIBRARIES ${FLENS_LIBRARY})
else (FLENS_FOUND)
  set (FLENS_INCLUDE_DIR)
  set (FLENS_LIBRARIES)
endif (FLENS_FOUND)
