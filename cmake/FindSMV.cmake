# - Find SMV
# Find the SMV headers and libraries
#
# SMV_INCLUDE_DIR - include path for headers
# SMV_LIBRARIES   - libraries to include when linking
# SMV_FOUND       - True if SMV is found

if (SMV_INCLUDE_DIR AND SMV_LIBRARIES)
  # already in cache, be silent
  set (SMV_FIND_QUIETLY TRUE)
endif (SMV_INCLUDE_DIR AND SMV_LIBRARIES)

# find the headers
find_path (SMV_INCLUDE_PATH
  SMV/Camera.h
  HINTS SMV_INC_DIR ENV SMV_INC_DIR
  )

find_library (SMV_LIBRARY
  NAMES SMV
  HINTS SMV_LIB_DIR ENV SMV_LIB_DIR
  )

# handle the QUIETLY and REQUIRED arguments and set
# SMV_FOUND to TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SMV "SMV not found. Set SMV_INC_DIR and SMV_LIB_DIR." SMV_INCLUDE_PATH SMV_LIBRARY)

if (SMV_FOUND)
  set (SMV_INCLUDE_DIR ${SMV_INCLUDE_PATH})
  set (SMV_LIBRARIES ${SMV_LIBRARY})
endif (SMV_FOUND)
