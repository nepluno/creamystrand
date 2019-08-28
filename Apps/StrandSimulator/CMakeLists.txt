include_directories (${CMAKE_CURRENT_SOURCE_DIR})

# Eigen 3 library is required
find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

# MKL is required
find_package( MKL )
if( MKL_FOUND)
    add_definitions(-DWITH_MKL)
    set ( LAPACK_INCLUDE_DIR ${MKL_INCLUDE_DIR} )
    set ( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
else( MKL_FOUND )
    set( ENV{BLA_VENDOR} ${SUGGESTED_BLAS_LAPACK} )
    find_package( LAPACK REQUIRED )
endif( MKL_FOUND )

find_package (TCLAP REQUIRED)
include_directories (${TCLAP_INCLUDE_DIR})

# Boost is required
# We need 1.58+ for the sort library
find_package( Boost 1.58.0 REQUIRED COMPONENTS serialization thread system filesystem )
include_directories( ${Boost_INCLUDE_DIRS} )

# OpenGL is required
find_package( OpenGL REQUIRED )
#include_directories( ${OPENGL_INCLUDE_DIR} )

# GLUT is required
find_package (GLUT REQUIRED glut)
if (GLUT_FOUND)
  include_directories (${GLUT_INCLUDE_DIR})
else (GLUT_FOUND)
  message (SEND_ERROR "Unable to locate GLUT")
endif (GLUT_FOUND)

# RapidXML library is required
find_package (RapidXML REQUIRED)
if (RAPIDXML_FOUND)
  include_directories (${RAPIDXML_INCLUDE_DIR})
else (RAPIDXML_FOUND)
  message (SEND_ERROR "Unable to locate RapidXML")
endif (RAPIDXML_FOUND)

if (USE_XCODE_PATCH AND CMAKE_GENERATOR STREQUAL Xcode)
  set (LLVM_ROOT_DIR "/usr/local/opt/llvm")
  find_package( LLVM REQUIRED )
  if( LLVM_FOUND )
    include_directories (${LLVM_LIBRARY_DIRS}/clang/${LLVM_VERSION_BASE_STRING}/include)
    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" )
  endif( LLVM_FOUND )
else()
  find_package( OpenMP )
  if( OPENMP_FOUND )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
  endif( OPENMP_FOUND )
endif()

# StrandSim test executable
AUX_SOURCE_DIRECTORY( Render RENDERDIR )
AUX_SOURCE_DIRECTORY( Problems PROBLEMDIR )

add_executable( StrandSimulatorApp StrandSimulatorApp.hh StrandSimulatorApp.cc  ProblemStepper.hh ProblemStepper.cc ${RENDERDIR} ${PROBLEMDIR} )
target_link_libraries( StrandSimulatorApp StrandSim bogus ${MKL_LIBRARY} ${LAPACK_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SERIALIZATION_LIBRARY} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY})

if (USE_XCODE_PATCH AND CMAKE_GENERATOR STREQUAL Xcode)
  set_target_properties(StrandSimulatorApp PROPERTIES XCODE_ATTRIBUTE_CC ${LLVM_ROOT_DIR}/bin/clang)
  set_target_properties(StrandSimulatorApp PROPERTIES XCODE_ATTRIBUTE_CXX ${LLVM_ROOT_DIR}/bin/clang)
  set_target_properties(StrandSimulatorApp PROPERTIES XCODE_ATTRIBUTE_COMPILER_INDEX_STORE_ENABLE "No")
endif()


# Creates a link to the assets directory in the binary output directory so the binary can access obj files and such
if (NOT WIN32)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/assets ${CMAKE_CURRENT_BINARY_DIR}/assets )
endif()





