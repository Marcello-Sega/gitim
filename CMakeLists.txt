CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11) 
PROJECT (ITIM CXX) 
INCLUDE (CheckIncludeFiles)
set_source_files_properties(gmx_itim.cpp PROPERTIES LANGUAGE CXX )
list( APPEND CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS}" )
if (NOT MYGMX_TREE OR NOT MYGMX_BUILD)
MESSAGE(FATAL_ERROR "MISSING MYGMX_TREE or MYGMX_BUILD options to specify the path of the gromacs source base directory (e.g. ~/gromacs-5.0/) and build directory (e.g. ~/gromacs-5.0/build/) . Pass them to cmake with -DMYGMX_TREE=..., etc. Note that if you compiled gromacs in the source directory, the two variable should point to the same directory.)")
ENDIF(NOT MYGMX_TREE OR NOT MYGMX_BUILD)

get_filename_component(MYGMX_BUILD ${MYGMX_BUILD}  REALPATH) 
get_filename_component(MYGMX_TREE ${MYGMX_TREE}  REALPATH) 

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB_RECURSE children "${curdir}/*")
  SET(dirlist "")
  FOREACH(child ${children}) 
    get_filename_component(child ${child} DIRECTORY)
    LIST(APPEND dirlist ${child})
  ENDFOREACH()
  LIST(REMOVE_DUPLICATES dirlist)
  SET(${result} ${dirlist})
ENDMACRO()


IF(GITIM_LR) 
 ADD_DEFINITIONS(-DGITIM_LR)
 MESSAGE(STATUS "Switching on suport for left/right interface (through the environment variables GITIM_LEFT and GITIM_RIGHT)" )
ENDIF(GITIM_LR)

# get_filename_component() is fundamental, as glob seems -so far- not to undersand '~'  

SET( GMX "${MYGMX_BUILD}" "${MYGMX_TREE}" )

# these are possible subdirectories where to find headers/libs, for each of the GMX dirs. This
# should take care of the tree differences between gmx v. 4 and v. 5, and minimize the changes in
# the code.
#SET (HEADER_SUBDIRS "include" "src" "lib" "src/gromacs" "src/external/boost/")

# let's combine everything

SET(MYPATHS ${MYGMX_TREE})
LIST(APPEND MYPATHS ${MYGMX_BUILD})
SUBDIRLIST(SUBDIRS ${MYGMX_TREE}/src)
LIST(APPEND MYPATHS ${SUBDIRS})
SUBDIRLIST(SUBDIRS ${MYGMX_TREE}/include)
LIST(APPEND MYPATHS ${SUBDIRS})
SUBDIRLIST(SUBDIRS ${MYGMX_BUILD})
LIST(APPEND MYPATHS ${SUBDIRS})
SUBDIRLIST(SUBDIRS ${MYGMX_BUILD}/lib)
LIST(APPEND MYPATHS ${SUBDIRS})

#FOREACH (DIR ${GMX})
#	FOREACH (SUBDIR ${HEADER_SUBDIRS})
#		LIST(APPEND MYPATHS ${DIR}${SUBDIR})
#	ENDFOREACH(SUBDIR)
#ENDFOREACH(DIR)


# let's first check which version are we using
# NOTE that gmx 4.6.x has another version.h file (without the info on the version...) in src/gmxlib/ 
#      we need to avoid that one by searching first for the include/ subdirs
FILE(GLOB_RECURSE versions "${MYGMX_BUILD}/*version.h")
IF(NOT versions)
	MESSAGE(FATAL_ERROR "Not able to find a version.h file")
ENDIF(NOT versions)
SET(GMX_VERSION NONE)
FOREACH(VERSION_H ${versions})
	MESSAGE(STATUS "checking version file ${VERSION_H}")
	FILE (READ ${VERSION_H} STRING_VERSION_H)
	IF("${STRING_VERSION_H}" MATCHES "#define GMX_VERSION" ) 
		STRING(REGEX REPLACE ".*#define GMX_VERSION ([0-9]+).*" "\\1" GMX_VERSION "${STRING_VERSION_H}")
		break()
	ENDIF()
ENDFOREACH(VERSION_H)

IF(${GMX_VERSION} MATCHES "NONE" )
        MESSAGE(FATAL_ERROR "Unable to determine gromacs version from file: ${VERSION_H}")
ENDIF()

IF(${UNIX}) 
# BUG: travis complaints about signal() are triggered only using GMX 5.1
IF(${GMX_VERSION} VERSION_LESS 50100) 
 ADD_DEFINITIONS(-DUNIX)
ENDIF()
ENDIF(${UNIX})


#### SET (MINGMXVERS 40601)
#### IF(${GMX_VERSION} VERSION_LESS ${MINGMXVERS})
####    MESSAGE(FATAL_ERROR "Version ${GMX_VERSION} < ${MINGMXVERS}")
#### ENDIF(${GMX_VERSION} VERSION_LESS ${MINGMXVERS})

MESSAGE(STATUS "============================================================")
MESSAGE(STATUS "GMX version: ${GMX_VERSION}    on    ${CMAKE_HOST_SYSTEM}   ")
ADD_DEFINITIONS(-DGMX_VERSION=${GMX_VERSION})
MESSAGE(STATUS "============================================================")

### version 4.6.x

IF(${GMX_VERSION} VERSION_LESS 50000) 
	# we need to check only for these two headers, the others will be in the same directories
	SET (HEADERLIST "gmx_statistics.h" "gmx_header_config.h" "version.h")
	SET (STATIC_LIBFILE libgmx.a libgmx_d.a)
	IF( CMAKE_HOST_SYSTEM MATCHES Darwin)
		SET(CMAKE_FIND_FRAMEWORK LAST)
		SET(DYNAMIC_LIBFILE libgmx.dylib libgmx_d.dylib)
	ELSE( CMAKE_HOST_SYSTEM MATCHES Darwin)
		SET (DYNAMIC_LIBFILE libgmx.so libgmx_d.so)
	ENDIF( CMAKE_HOST_SYSTEM MATCHES Darwin)
ELSE(${GMX_VERSION} VERSION_LESS 50000) 
### version 5.x.x
	SET (HEADERLIST "config.h" "version.h" "smalloc.h" "gstat.h" "vec.h" "xvgr.h" "pbc.h" "copyrite.h" "futil.h" "tpxio.h" "gmx_ana.h" "nbsearch.h" "gromacs/commandline/cmdlinemodulemanager.h" "utility.h" "vec.h" "shared_ptr.hpp" "baseversion-gen.c")
	IF(${GMX_VERSION} VERSION_GREATER 50000) 
		SET (HEADERLIST ${HEADERLIST} "index.h" )
	ENDIF(${GMX_VERSION} VERSION_GREATER 50000)
	SET (STATIC_LIBFILE libgromacs.a libgromacs_d.a)
	IF( CMAKE_HOST_SYSTEM MATCHES Darwin)
		SET(CMAKE_FIND_FRAMEWORK LAST)
		SET (DYNAMIC_LIBFILE libgromacs.dylib libgromacs_d.dylib)
	ELSE( CMAKE_HOST_SYSTEM MATCHES Darwin)
		SET (DYNAMIC_LIBFILE libgromacs.so libgromacs_d.so)
	ENDIF( CMAKE_HOST_SYSTEM MATCHES Darwin)

ENDIF(${GMX_VERSION} VERSION_LESS 50000) 

### find libraries and determine if the version we need is the single or double precision
SET(GROMACS_LIBGMX_STATIC NOTFOUND)
SET(GROMACS_LIBGMX_DYNAMIC NOTFOUND)

### Static libraries
## search for double precision libraries
#FIND_PATH(GROMACS_LIBGMX_STATIC NAMES ${STATIC_LIBFILE} PATHS ${MYPATHS})
MESSAGE(STATUS "searching for ${STATIC_LIBFILE}")
FOREACH(SUBDIR ${MYPATHS})
	FOREACH(LIB ${STATIC_LIBFILE})
		FILE(GLOB GROMACS_LIBGMX_STATIC ${SUBDIR}/${LIB})
		IF(GROMACS_LIBGMX_STATIC)
			LIST(GET GROMACS_LIBGMX_STATIC 0 GROMACS_LIBGMX_STATIC)
			MESSAGE(STATUS "found ${GROMACS_LIBGMX_STATIC}")
			BREAK()
		ENDIF(GROMACS_LIBGMX_STATIC)
	ENDFOREACH(LIB ${STATIC_LIBFILE})
	IF(GROMACS_LIBGMX_STATIC)
			BREAK()
	ENDIF(GROMACS_LIBGMX_STATIC)
ENDFOREACH(SUBDIR ${MYPATHS})

MESSAGE(STATUS "searching for ${DYNAMIC_LIBFILE}")
FOREACH(SUBDIR ${MYPATHS})
	FOREACH(LIB ${DYNAMIC_LIBFILE})
		FILE(GLOB GROMACS_LIBGMX_DYNAMIC ${SUBDIR}/${LIB})
		IF(GROMACS_LIBGMX_DYNAMIC)
			LIST(GET GROMACS_LIBGMX_DYNAMIC 0 GROMACS_LIBGMX_DYNAMIC)
			MESSAGE(STATUS "found ${GROMACS_LIBGMX_DYNAMIC}")
			BREAK()
		ENDIF(GROMACS_LIBGMX_DYNAMIC)
	ENDFOREACH(LIB ${STATIC_LIBFILE})
	IF(GROMACS_LIBGMX_DYNAMIC)
			BREAK()
	ENDIF(GROMACS_LIBGMX_DYNAMIC)
ENDFOREACH(SUBDIR ${MYPATHS})



IF(NOT GROMACS_LIBGMX_STATIC  AND NOT GROMACS_LIBGMX_DYNAMIC )
       MESSAGE(FATAL_ERROR "Unable to find gromacs libraries")
ENDIF(NOT GROMACS_LIBGMX_STATIC  AND NOT GROMACS_LIBGMX_DYNAMIC )


IF(GROMACS_LIBGMX_STATIC)
	SET(GROMACS_LIBGMX "${GROMACS_LIBGMX_STATIC}")
ENDIF(GROMACS_LIBGMX_STATIC)
IF(GROMACS_LIBGMX_DYNAMIC)
	SET(GROMACS_LIBGMX "${GROMACS_LIBGMX_DYNAMIC}")
ENDIF(GROMACS_LIBGMX_DYNAMIC)
MESSAGE(STATUS "Using gromacs libraries: ${GROMACS_LIBGMX}")

        
SET (SUFFIXES gmxana gmxlib math fileio pbcutil legacyheaders selection utility pbcutil boost)

FOREACH(HEADER ${HEADERLIST})
	# looks like we have to reset TMPDIR every time to the failure value.
	SET(TMPDIR "TMPDIR-NOTFOUND")
	FIND_PATH(TMPDIR ${HEADER} PATHS ${MYPATHS} NO_DEFAULT_PATH PATH_SUFFIXES ${SUFFIXES})
	IF(NOT TMPDIR)
   		MESSAGE(FATAL_ERROR "Unable to find ${HEADER} in ${MYPATHS} :: ${TMPDIR}")
	ENDIF(NOT TMPDIR)
	FOREACH(SUFFIX ${SUFFIXES})
		STRING(FIND ${TMPDIR} ${SUFFIX} FOUND )
		IF(${FOUND} GREATER -1 )
		INCLUDE_DIRECTORIES(${TMPDIR}/../)
		ENDIF()
	ENDFOREACH()

	INCLUDE_DIRECTORIES(${TMPDIR})
ENDFOREACH(HEADER)
IF(NOT GROMACS_LIBGMX)
   MESSAGE(FATAL_ERROR "Unable to find libraries")
ENDIF(NOT GROMACS_LIBGMX)

## check source for gromacs helpcontext interface
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
SET(CMAKE_REQUIRED_INCLUDES "${inc_dirs}")
SET(CMAKE_REQUIRED_LIBRARIES "${GROMACS_LIBGMX}")
INCLUDE(CheckCXXSourceCompiles) 

IF(NOT ${GMX_VERSION} VERSION_LESS 50000) 
unset(HAVE_NEW_GMX_CONTEXT_INTERFACE CACHE)
unset(HAVE_OLD_GMX_CONTEXT_INTERFACE CACHE)
CHECK_CXX_SOURCE_COMPILES(
"
#include \"gromacs/commandline/cmdlinemodulemanager.h\"
#include \"gromacs/commandline/cmdlinehelpcontext.h\"
#include \"file.h\"
int
main()
{
            gmx::CommandLineHelpContext context(&gmx::File::standardError(),
			gmx::eHelpOutputFormat_Console, NULL,\"itim\");
}
" HAVE_NEW_GMX_CONTEXT_INTERFACE)


CHECK_CXX_SOURCE_COMPILES(
"
#include \"gromacs/commandline/cmdlinemodulemanager.h\"
#include \"gromacs/commandline/cmdlinehelpcontext.h\"
#include \"file.h\"
int
main()
{
            gmx::CommandLineHelpContext context(&gmx::File::standardError(),
			gmx::eHelpOutputFormat_Console, NULL);
}
" HAVE_OLD_GMX_CONTEXT_INTERFACE)


IF(NOT HAVE_NEW_GMX_CONTEXT_INTERFACE)
  IF(NOT HAVE_OLD_GMX_CONTEXT_INTERFACE)
    MESSAGE(FATAL_ERROR "No working interface for help found")
  ENDIF(NOT HAVE_OLD_GMX_CONTEXT_INTERFACE)
ENDIF(NOT HAVE_NEW_GMX_CONTEXT_INTERFACE)


IF(HAVE_NEW_GMX_CONTEXT_INTERFACE)
 ADD_DEFINITIONS(-DNEW_CONTEXT_INTERFACE)
ENDIF(HAVE_NEW_GMX_CONTEXT_INTERFACE)
ENDIF(NOT ${GMX_VERSION} VERSION_LESS 50000) 
	

MESSAGE(STATUS "")
MESSAGE(STATUS "=========================================================================")
MESSAGE(STATUS "Using libraries: ${GROMACS_LIBGMX} ")
MESSAGE(STATUS "=========================================================================")
MESSAGE(STATUS "")


# Let's tell cmake to check parameters for the examples
ADD_DEFINITIONS(-Wall  -pedantic  -Wno-unused -Wuninitialized )
IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
ADD_DEFINITIONS( -O0 )
ELSE (CMAKE_BUILD_TYPE STREQUAL "Debug")
ADD_DEFINITIONS( -O3 -ffast-math )
ENDIF(CMAKE_BUILD_TYPE STREQUAL "Debug")


ADD_EXECUTABLE(itim gmx_itim.cpp)
TARGET_LINK_LIBRARIES(itim pthread  ${GROMACS_LIBGMX} m dl )

IF(${GMX_VERSION} VERSION_LESS 50000)
# todo
ELSE(${GMX_VERSION} VERSION_LESS 50000)
	FILE (WRITE ${CMAKE_SOURCE_DIR}/examples/ccl4-h2o/include.sh "GROMPP='${MYGMX_BUILD}/bin/gmx grompp'\nMAKENDX='${MYGMX_BUILD}/bin/gmx make_ndx'\nGDENS=../../itim\nTOPOL=topol-5.0.top\nexport GMXLIB=${MYGMX_TREE}/share/top/\n" )
ENDIF(${GMX_VERSION} VERSION_LESS 50000)

add_custom_command(TARGET  itim PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/examples/ccl4-h2o $<TARGET_FILE_DIR:itim>/examples/ccl4-h2o)


