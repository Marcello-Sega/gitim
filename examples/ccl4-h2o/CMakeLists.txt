CMAKE_MINIMUM_REQUIRED(VERSION 2.6) 
PROJECT (GITIM C) 
MESSAGE(STATUS "Entering examples/ccl4-h2o")
FILE(WRITE included.sh "export GMXLIB=${GROMACS_TOP_DIR}\nexport GROMPP=${GROMACS_KERNEL_DIR}/grompp \nexport GDENS=../../gitim\n")


