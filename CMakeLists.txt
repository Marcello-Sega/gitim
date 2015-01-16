CMAKE_MINIMUM_REQUIRED(VERSION 2.6) 
PROJECT (G_DENSITY CXX) 
INCLUDE (CheckIncludeFiles)
set_source_files_properties(g_density.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(gmx_density.c PROPERTIES LANGUAGE CXX )
list( APPEND CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS}" )
SET(ENVSRCDIR $ENV{MYGMX_SRC})
SET(ENVBLDDIR $ENV{MYGMX_BUILD})
if(NOT ENVSRCDIR)
	MESSAGE(STATUS $ENV{MYGMX_SRC})
	MESSAGE(FATAL_ERROR "\nERROR: set the variables MYGMX_SRC  and MYGMX_BUILD to specify a particular path for gromacs source and build directories" )
endif(NOT ENVSRCDIR)
#IF(VIRIAL_EXTENSION)
#	 add_definitions(-DVIRIAL_EXTENSION)
#ENDIF(VIRIAL_EXTENSION)
# these are possible locatsions of gromacs' distribution source / headers / libraries
SET( GMX
"$ENV{MYGMX_BUILD}/"
"$ENV{MYGMX_SRC}/"
"/usr/include/"
"/usr/lib/"
"/usr/local/include/"
"/usr/local/lib/"
)
# these are possible subdirectories where to find headers/libs, for each of the GMX dirs. This
# should take care of the tree differences between gmx v. 4 and v. 5, and minimize the changes in
# the code.
SET (HEADER_SUBDIRS "include" "src/" "src/gmxlib"  "src/gromacs/" "src/gromacs/utility" "src/gromacs/gmxana/" "src/gromacs/math/" "src/gromacs/fileio/" "src/gromacs/math/" "src/gromacs/pbcutil/" "src/gromacs/legacyheaders/" "src/gromacs/selection/" "lib/" "gmxlib"  "gromacs/" "gromacs/utility" "gromacs/gmxana/" "gromacs/math/" "gromacs/fileio/" "gromacs/math/" "gromacs/pbcutil/" "gromacs/legacyheaders/" "gromacs/selection/" )

# let's combine everything
SET(MYPATHS    ${GMX})
FOREACH (DIR ${GMX})
	FOREACH (SUBDIR ${HEADER_SUBDIRS})
		LIST(APPEND MYPATHS "${DIR}/${SUBDIR}/")
	ENDFOREACH(SUBDIR)
ENDFOREACH(DIR)

# let's first check which version are we using
# NOTE that gmx 4.6.x has another version.h file (without the info on the version...) in src/gmxlib/ 
#      we need to avoid that one by searching first for the include/ subdirs
FIND_FILE(VERSION_H version.h PATHS ${MYPATHS})
IF(NOT VERSION_H)
	MESSAGE(FATAL_ERROR "Not able to find a version.h file in ${MYPATHS}")
ENDIF(NOT VERSION_H)

MESSAGE(STATUS "Using version file ${VERSION_H}")
FILE (STRINGS ${VERSION_H} STRING_VERSION_H)
IF(NOT  "${STRING_VERSION_H}" MATCHES "#define GMX_VERSION" ) 
	MESSAGE(FATAL_ERROR "Unable to determine gromacs version from file: ${VERSION_H}")
ENDIF(NOT  "${STRING_VERSION_H}" MATCHES "#define GMX_VERSION" ) 

STRING(REGEX REPLACE ".*#define GMX_VERSION ([0-9]+).*" "\\1" GMX_VERSION "${STRING_VERSION_H}")
SET (MINGMXVERS 40601)
IF(${GMX_VERSION} VERSION_LESS ${MINGMXVERS})
   MESSAGE(FATAL_ERROR "Version ${GMX_VERSION} < ${MINGMXVERS}")
ENDIF(${GMX_VERSION} VERSION_LESS ${MINGMXVERS})

MESSAGE(STATUS "============================================================")
MESSAGE(STATUS "GMX version: ${GMX_VERSION}    on    ${CMAKE_HOST_SYSTEM}   ")
ADD_DEFINITIONS(-DGMX_VERSION=${GMX_VERSION})
MESSAGE(STATUS "============================================================")

### version 4.6.x

IF(${GMX_VERSION} VERSION_LESS 50000) 
	# we need to check only for these two headers, the others will be in the same directories
	SET (HEADERLIST "gmx_statistics.h" "gmx_header_config.h" "version.h")
	SET (STATIC_LIBFILE libgmx.a)
	IF( CMAKE_HOST_SYSTEM MATCHES Darwin)
		SET(DYNAMIC_LIBFILE libgmx.dylib)
	ELSE( CMAKE_HOST_SYSTEM MATCHES Darwin)
		SET (DYNAMIC_LIBFILE libgmx.so)
	ENDIF( CMAKE_HOST_SYSTEM MATCHES Darwin)
ELSE(${GMX_VERSION} VERSION_LESS 50000) 

### version 5.x.x
	SET (HEADERLIST "version.h" "string.h" "smalloc.h" "gstat.h" "vec.h" "xvgr.h" "pbc.h" "copyrite.h" "futil.h" "tpxio.h" "gmx_ana.h" "nbsearch.h" "gmxpre.h")
	SET (STATIC_LIBFILE libgromacs.a)
	IF( CMAKE_HOST_SYSTEM MATCHES Darwin)
		SET (DYNAMIC_LIBFILE libgromacs.dylib)
	ELSE( CMAKE_HOST_SYSTEM MATCHES Darwin)
		SET (DYNAMIC_LIBFILE libgromacs.so)
	ENDIF( CMAKE_HOST_SYSTEM MATCHES Darwin)
ENDIF(${GMX_VERSION} VERSION_LESS 50000) 

### endof version-dependent part

FIND_PATH(GROMACS_LIBGMX_STATIC NAMES ${STATIC_LIBFILE} PATHS ${MYPATHS})
FIND_PATH(GROMACS_LIBGMX_DYNAMIC NAMES ${DYNAMIC_LIBFILE} PATHS ${MYPATHS})

IF(NOT GROMACS_LIBGMX_STATIC AND NOT GROMACS_LIBGMX_DYNAMIC)
       MESSAGE(FATAL_ERROR "Unable to find gromacs libraries")
ENDIF(NOT GROMACS_LIBGMX_STATIC  AND NOT GROMACS_LIBGMX_DYNAMIC) 

IF(GROMACS_LIBGMX_STATIC)
	SET(GROMACS_LIBGMX "${GROMACS_LIBGMX_STATIC}/${STATIC_LIBFILE}")
ENDIF(GROMACS_LIBGMX_STATIC)
IF(GROMACS_LIBGMX_DYNAMIC)
	SET(GROMACS_LIBGMX "${GROMACS_LIBGMX_DYNAMIC}/${DYNAMIC_LIBFILE}")
ENDIF(GROMACS_LIBGMX_DYNAMIC)



FOREACH(HEADER ${HEADERLIST})
	# looks like we have to reset TMPDIR every time to the failure value.
	SET(TMPDIR "TMPDIR-NOTFOUND")
	FIND_PATH(TMPDIR ${HEADER} PATHS ${MYPATHS})
	IF(NOT TMPDIR)
   		MESSAGE(FATAL_ERROR "Unable to find ${HEADER} in ${MYPATHS} :: ${TMPDIR}")
	ENDIF(NOT TMPDIR)
	MESSAGE(STATUS "found ${HEADER} in ${TMPDIR}")
	INCLUDE_DIRECTORIES(${TMPDIR})
ENDFOREACH(HEADER)
IF(NOT GROMACS_LIBGMX)
   MESSAGE(FATAL_ERROR "Unable to find libraries")
ENDIF(NOT GROMACS_LIBGMX)
	
#-I/home/marcello/gromacs.gerrit/src/gromacs/selection -I /home/marcello/gromacs.gerrit/src/gromacs/legacyheaders/ -I/home/marcello/gromacs.gerrit/src/gromacs/pbcutil/ -I/home/marcello/gromacs.gerrit/src/gromacs/fileio/  -I/home/marcello/gromacs.gerrit/src/gromacs/math/  -I/home/marcello/gromacs.gerrit/src -I/home/marcello/gromacs.gerrit/src/gromacs/utility/  -I/home/marcello/gromacs.gerrit/build/  -I/home/marcello/gromacs.gerrit/src/gromacs/gmxana/





##get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
##MESSAGE(STATUS "Found gromacs headers in ${dirs}")

MESSAGE(STATUS "")
MESSAGE(STATUS "==================================================")
MESSAGE(STATUS "Using libraries in: ${GROMACS_LIBGMX} ")
MESSAGE(STATUS "Using headers in: ${TMPDIR}")
MESSAGE(STATUS "==================================================")
MESSAGE(STATUS "")



# Let's tell cmake to check parameters for the examples
ADD_SUBDIRECTORY(examples/ccl4-h2o)
GET_DIRECTORY_PROPERTY(GROMACS_TOP_DIR DIRECTORY examples/ccl4-h2o DEFINITION GROMACS_TOP_DIR)
GET_DIRECTORY_PROPERTY(GROMACS_KERNEL_DIR DIRECTORY examples/ccl4-h2o DEFINITION GROMACS_KERNEL_DIR)

ADD_DEFINITIONS( -O3 -Wall  -pedantic  -Wno-unused -Wuninitialized  -ffast-math )
ADD_EXECUTABLE(g_density g_density.c gmx_density.c)
MESSAGE(STATUS "GMX LIB:  ${GROMACS_LIBGMX}")
TARGET_LINK_LIBRARIES(g_density pthread  ${GROMACS_LIBGMX} m dl )



