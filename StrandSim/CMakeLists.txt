include( CMakeSourceFiles.txt )

# Eigen 3 library is required
find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

# Approximate Nearest Neighbors is required
find_package( ANN REQUIRED )
include_directories( ${ANN_INCLUDE_DIR} )

# MKL is required
set (MKL_MULTI_THREADED ON)
find_package( MKL )
if( MKL_FOUND)
    add_definitions(-DWITH_MKL)
    set ( LAPACK_INCLUDE_DIR ${MKL_INCLUDE_DIR} )
    set ( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
else( MKL_FOUND )
    set( ENV{BLA_VENDOR} ${SUGGESTED_BLAS_LAPACK} )
    find_package( LAPACK REQUIRED )
endif( MKL_FOUND )

include_directories( ${LAPACK_INCLUDE_DIR} )

# Boost is required
# We need 1.58+ for the sort library
find_package( Boost 1.58.0 REQUIRED COMPONENTS serialization thread system filesystem )
include_directories( ${Boost_INCLUDE_DIRS} )

if (USE_XCODE_PATCH AND CMAKE_GENERATOR STREQUAL Xcode)
	set (LLVM_ROOT_DIR "/usr/local/opt/llvm")
	find_package( LLVM REQUIRED )
	if( LLVM_FOUND )
		include_directories (${LLVM_LIBRARY_DIRS}/clang/${LLVM_VERSION_BASE_STRING}/include)
		set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" )
    add_definitions(-D_OPENMP)
	endif( LLVM_FOUND )
else()
	find_package( OpenMP )
	if( OPENMP_FOUND )
		set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
  else()
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS}")
	endif( OPENMP_FOUND )
endif()

# IGL library is required
find_package (LIBIGL REQUIRED)
if (LIBIGL_FOUND)
  include_directories (${LIBIGL_INCLUDE_DIR})
else (LIBIGL_FOUND)
  message (SEND_ERROR "Unable to locate IGL")
endif (LIBIGL_FOUND)

# GLUT is required
find_package (GLUT REQUIRED glut)
if (GLUT_FOUND)
  include_directories (${GLUT_INCLUDE_DIR})
else (GLUT_FOUND)
  message (SEND_ERROR "Unable to locate GLUT")
endif (GLUT_FOUND)

# StrandSim library
add_library( StrandSim ${Headers} ${Sources} )

if (USE_XCODE_PATCH AND CMAKE_GENERATOR STREQUAL Xcode)
  set_target_properties(StrandSim PROPERTIES XCODE_ATTRIBUTE_CC ${LLVM_ROOT_DIR}/bin/clang)
  set_target_properties(StrandSim PROPERTIES XCODE_ATTRIBUTE_CXX ${LLVM_ROOT_DIR}/bin/clang)
  set_target_properties(StrandSim PROPERTIES XCODE_ATTRIBUTE_COMPILER_INDEX_STORE_ENABLE "No")
endif()
