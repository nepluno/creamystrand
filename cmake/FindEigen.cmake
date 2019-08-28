# - Find Eigen
# Find the Eigen headers
#
# EIGEN_INCLUDE_DIR - where to find the Eigen headers
# EIGEN_FOUND       - True if Eigen is found

if (EIGEN_INCLUDE_DIR)
  # already in cache, be silent
  set (EIGEN_FIND_QUIETLY TRUE)
endif (EIGEN_INCLUDE_DIR)

# find the headers
find_path (EIGEN_INCLUDE_PATH
  Eigen/Core
  HINTS Eigen_INC_DIR ENV Eigen_INC_DIR
  DOC "Eigen include directory (http://eigen.tuxfamily.org/index.php?title=Main_Page)"
  )

# handle the QUIETLY and REQUIRED arguments and set EIGEN_FOUND to
# TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Eigen "Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page) could not be found. Set Eigen_INC_DIR to point to the headers." EIGEN_INCLUDE_PATH)

if (EIGEN_FOUND)
  set (EIGEN_INCLUDE_DIR ${EIGEN_INCLUDE_PATH})
endif (EIGEN_FOUND)
