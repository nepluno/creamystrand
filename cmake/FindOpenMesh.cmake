# - Find OpenMesh
# Find the OpenMesh includes and library
#
# OPENMESH_INCLUDE_DIR - where to find the OpenMesh headers
# OPENMESH_LIBRARIES   - list of libraries when using OpenMesh
# OPENMESH_FOUND       - True if OpenMesh is found

if (OPENMESH_INCLUDE_DIR AND OPENMESH_LIBRARIES)
  # already in cache, be silent
  set (OPENMESH_FIND_QUIETLY TRUE)
endif (OPENMESH_INCLUDE_DIR AND OPENMESH_LIBRARIES)

# find the headers
find_path (OPENMESH_INCLUDE_PATH
  OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh
  HINTS OpenMesh_INC_DIR $ENV{OpenMesh_INC_DIR}
  )

# find the libraries
set (OPENMESH_NAMES OpenMeshCore)
find_library (OPENMESH_LIBRARY NAMES
  NAMES ${OPENMESH_NAMES}
  HINTS OpenMesh_LIB_DIR $ENV{OpenMesh_LIB_DIR}
  )

# handle the QUIETLY and REQUIRED arguments and set OPENMESH_FOUND to
# TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (OpenMesh "OpenMesh (http://openmesh.org/) not found. Set OpenMesh_INC_DIR and OpenMesh_LIB_DIR." OPENMESH_INCLUDE_PATH OPENMESH_LIBRARY)

if (OPENMESH_FOUND)
  set (OPENMESH_INCLUDE_DIR ${OPENMESH_INCLUDE_PATH})
  set (OPENMESH_LIBRARIES ${OPENMESH_LIBRARY})
else (OPENMESH_FOUND)
  set (OPENMESH_INCLUDE_DIR)
  set (OPENMESH_LIBRARIES)
endif (OPENMESH_FOUND)
