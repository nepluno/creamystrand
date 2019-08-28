# - Find OpenCurrent
# Find the OpenCurrent headers and libraries
#
# OCU_INCLUDE_DIR - include path for headers
# OCU_LIBRARIES   - libraries to include when linking
# OCU_FOUND       - True if MAYA is found

if (OCU_INCLUDE_DIR AND OCU_LIBRARIES)
  # already in cache, be silent
  set (OCU_FIND_QUIETLY TRUE)
endif (OCU_INCLUDE_DIR AND OCU_LIBRARIES)


##############################################################
##
## Search for the header location. There are hundreds of
## headers so just search for one we know must be there and 
## assume it is a valid install.
##
##############################################################
find_path (OCU_INCLUDE_PATH
  /ocustorage/grid3d.h
  HINTS ENV OCU_LOCATION
  )

##############################################################
##
## Search for all the libraries to make sure they can be found
##
##############################################################
set (OCU_LIBRARIES_FOUND)
set (OCU_LIBRARIES_MISSING)
set (OCU_LIBS ocuequation ocustorage ocuutil)
foreach(OCULIB ${OCU_LIBS})
	       set (OCU_SEARCH_LIB "OCU_SEARCH_LIB-NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
	       find_library (OCU_SEARCH_LIB ${OCULIB} PATHS HINTS $ENV{OCU_LOCATION}/sm13-rel/${OCULIB})
	       if(OCU_SEARCH_LIB)
		    list (APPEND OCU_LIBRARIES_FOUND ${OCU_SEARCH_LIB})
	       else (OCU_SEARCH_LIB)
		    list (APPEND OCU_LIBRARIES_MISSING ${OCULIB})
		    message (SEND_ERROR "Unable to find Opencurrent library ${OCULIB}")
	       endif (OCU_SEARCH_LIB)
endforeach(OCULIB)

set (OCU_LIBRARY ${OCU_LIBRARIES_FOUND} CACHE STRING "Opencurrent libraries" FORCE)


# handle the QUIETLY and REQUIRED arguments and set
# OCU_FOUND to TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
#find_package_handle_standard_args (OCU "OCU not found. Set OCU_LOCATION " OCU_INCLUDE_PATH OCU_LIBRARIES)
find_package_handle_standard_args (OCU "OCU not found. Set OCU_LOCATION " OCU_INCLUDE_PATH OCU_LIBRARY)

if (OCU_FOUND)
  set (OCU_INCLUDE_DIR ${OCU_INCLUDE_PATH}/src)
  set (OCU_LIBRARIES ${OCU_LIBRARY})
else (OCU_FOUND)
  set (OCU_INCLUDE_DIR)
  set (OCU_LIBRARIES)
endif (OCU_FOUND)
