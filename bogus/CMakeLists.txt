
execute_process( COMMAND ../updateCMakeSources.sh WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
file(STRINGS CMakeSources.txt SRC_LIST)

find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

find_package( Boost COMPONENTS serialization thread system )

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()

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
	endif( OPENMP_FOUND )
endif()

add_definitions( -DBOGUS_INSTANTIATE_3D_SOC)

if( WITH_DYNAMIC )
        add_definitions( -DBOGUS_INSTANTIATE_DYNAMIC_SOC )
endif()

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -S --verbose-asm" )

add_library( bogus ${SRC_LIST})

if (USE_XCODE_PATCH AND CMAKE_GENERATOR STREQUAL Xcode)
  set_target_properties(bogus PROPERTIES XCODE_ATTRIBUTE_CC ${LLVM_ROOT_DIR}/bin/clang)
  set_target_properties(bogus PROPERTIES XCODE_ATTRIBUTE_CXX ${LLVM_ROOT_DIR}/bin/clang)
  set_target_properties(bogus PROPERTIES XCODE_ATTRIBUTE_COMPILER_INDEX_STORE_ENABLE "No")
endif()

FOREACH(item ${SRC_LIST})
  STRING(REGEX MATCH "^Interfaces/.+\\.hpp$" item ${item})
  IF(item)
	LIST(APPEND INTERFACES ${item})
  ENDIF(item)
ENDFOREACH(item ${SRC_LIST})


