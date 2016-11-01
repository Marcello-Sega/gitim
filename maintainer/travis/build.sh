#!/usr/bin/env bash
# Copyright (C) 2015 Marcello Sega
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.


# HELPER FUNCTIONS

# output value of env variables
function outp {
    for p in $*; do
        echo "  $p=${!p}"
    done
}

# start a block
function start {
    echo "=================================================="
    echo "START $1"
    echo "=================================================="
}

# end a block
function end {
    echo "=================================================="
    echo "END $1"
    echo "=================================================="
}

# execute and output a command
function cmd {
    echo ">$1"
    eval $1
}


# handle environment variables
[ -z "$insource" ] && insource="true"
[ -z "$srcdir" ] && srcdir=`pwd`
#[ -z "$configure_params" ] && configure_params=""
#[ -z "$configure_vars" ] && configure_vars=""
#[ -z "$with_mpi" ] && with_mpi="true"
#[ -z "$with_fftw" ] && with_fftw="true"
#[ -z "$with_tcl" ] && with_tcl="true"
#[ -z "$with_python_interface" ] && with_python_interface="true"
#[ -z "$myconfig" ] && myconfig="default"
#! $with_mpi && check_procs=1
#[ -z "$check_procs" ] && check_procs=4
#[ -z "$make_check" ] && make_check="true"

if $insource; then
    builddir=$srcdir
elif [ -z "$builddir" ]; then
    builddir=$srcdir/build
fi

outp insource srcdir builddir 

if ! $insource; then
    if [ ! -d $builddir ]; then
        echo "Creating $builddir..."
        mkdir -p $builddir
    fi
fi

# BOOTSTRAP
start "GROMACS DOWNLOAD"
cmd "curl -O ftp://ftp.gromacs.org/pub/gromacs/gromacs-${GMX}.tar.gz"
end "GROMACS DOWNLOAD"

if ! $insource ; then
    cd $builddir
fi

start "GROMACS BUILD"
cmd "tar -xzf gromacs-${GMX}.tar.gz"
cmd "mv gromacs-${GMX} gromacs"
cmd "mkdir gromacs/build"
(
 cmd "cd gromacs/build" || exit $?
 cmd "cmake .. -DGMX_FFT_LIBRARY=fftpack -DGMX_GPU=OFF -DGMX_SIMD=SSE4.1" || exit $? 
 cmd "make" || exit $? 
) 
end "GROMACS BUILD"
