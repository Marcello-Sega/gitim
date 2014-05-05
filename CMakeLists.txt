CMAKE_MINIMUM_REQUIRED(VERSION 2.6) 
PROJECT (G_DENSITY C) 
FILE(GLOB GMX 
"$ENV{HOME}/gromacs-*/include" 
"/usr/local/gromacs-*/include"
"/usr/gromacs-*/include" 
"/opt/gromacs-*/include")

SET(MYPATHS      /usr/include 
                $ENV{HOME}/include 
                /usr/local/include
                ${GMX})

MESSAGE(STATUS "searching in ${MYPATHS}")
FIND_PATH(GROMACS_INC_DIR NAMES gmx_statistics.h PATHS ${MYPATHS})

IF(NOT GROMACS_INC_DIR)
   MESSAGE(FATAL_ERROR "Unable to find gromacs include files")
ELSE(NOT GROMACS_INC_DIR)
   MESSAGE(STATUS "Found gromacs headers in ${GROMACS_INC_DIR}")
   STRING(REGEX REPLACE "/include" "/share/top" GROMACS_TOP_DIR "${GROMACS_INC_DIR}")
   STRING(REGEX REPLACE "/include" "/src/kernel" GROMACS_KERNEL_DIR "${GROMACS_INC_DIR}")
   STRING(REGEX REPLACE "/include" "/src/gmxlib/.libs/" GROMACS_LIB_DIR "${GROMACS_INC_DIR}")
   FIND_PATH(GROMACS_LIBGMX NAMES libgmx.a PATHS ${GROMACS_LIB_DIR})
   IF(NOT GROMACS_LIBGMX)
       MESSAGE(STATUS "Unable to find libgmx.a in ${GROMACS_LIBGMX}")
   ELSE(NOT GROMACS_LIBGMX)
       MESSAGE(STATUS "Found libgmx.a in ${GROMACS_LIBGMX}")
       SET(GROMACS_LIBGMX "${GROMACS_LIBGMX}/libgmx.a")
   ENDIF(NOT GROMACS_LIBGMX) 
ENDIF(NOT GROMACS_INC_DIR) 




INCLUDE_DIRECTORIES(${GROMACS_INC_DIR})
# Let's tell cmake to check parameters for the examples
ADD_SUBDIRECTORY(examples/ccl4-h2o)
GET_DIRECTORY_PROPERTY(GROMACS_TOP_DIR DIRECTORY examples/ccl4-h2o DEFINITION GROMACS_TOP_DIR)
GET_DIRECTORY_PROPERTY(GROMACS_KERNEL_DIR DIRECTORY examples/ccl4-h2o DEFINITION GROMACS_KERNEL_DIR)

ADD_DEFINITIONS( -O3 -Wall  -pedantic -std=c99  -Wno-unused -Wuninitialized  -ffast-math )
ADD_DEFINITIONS( -DALPHA ) 
ADD_EXECUTABLE(g_density g_density.c gmx_density.c)
TARGET_LINK_LIBRARIES(g_density pthread  ${GROMACS_LIBGMX} m dl )



