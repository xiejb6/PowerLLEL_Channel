#!/bin/bash

# Description: This script shows how to build PowerLLEL and run a demo on Tianhe-2A.
#              Make necessary modifications to port to your machine.

# 1. Load environments.
# 1) By package management tool.
#   Spack is a handy package management tool designed for large supercomputing centers
#   (https://github.com/spack/spack). Another similar tool is Environment Modules 
#   (https://modules.readthedocs.io/).
# 2) Manually.
#   Set environment variables ($PATH, $LD_LIBRARY_PATH, $LIBRARY_PATH, $CPATH, $PKG_CONFIG_PATH, etc.)
#   directly using the following commands:
#     export PATH=${PATH_TO_ADD}:$PATH
#     export LD_LIBRARY_PATH=${LD_LIBRARY_PATH_TO_ADD}:$LD_LIBRARY_PATH
#     ...

set -e

source ~/fgn/env/cn/intel-18.0.0.sh

# LOAD HDF5
export PATH=$HOME/fgn/software/hdf5/1.10.4/build/install/bin:$PATH
export LD_LIBRARY_PATH=$HOME/fgn/software/hdf5/1.10.4/build/install/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$HOME/fgn/software/hdf5/1.10.4/build/install/lib:$LIBRARY_PATH
export INCLUDE=$HOME/fgn/software/hdf5/1.10.4/build/install/include:$INCLUDE

# LOAD FFTW
export PKG_CONFIG_PATH=$HOME/fgn/software/fftw/3.3.8/build/install/lib/pkgconfig:$PKG_CONFIG_PATH

# LOAD GPTL
export LD_LIBRARY_PATH=$HOME/fgn/software/gptl/8.1.1/build/install/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$HOME/fgn/software/gptl/8.1.1/build/install/lib/pkgconfig:$PKG_CONFIG_PATH

# 2. Modify CMake options shown below as needed

rm -rf build
cmake \
    -DCMAKE_INSTALL_PREFIX=build/install \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_C_COMPILER=icc \
    -DMPI_C_COMPILER=mpicc \
    -DCMAKE_C_FLAGS="-O3 -xHost -ipo" \
    -DCMAKE_Fortran_COMPILER=mpifort \
    -DMPI_Fortran_COMPILER=mpifort \
    -DCMAKE_Fortran_FLAGS="-O3 -xHost -ipo" \
    -DMPI_Fortran_HAVE_F08_MODULE=false \
    `# IMPORTANT options` \
    -DUSE_GPTL=ON \
    -DUSE_C=OFF \
    -DUSE_OMP=OFF \
    -DUSE_PDD=OFF \
    -DUSE_MKL=OFF \
    -DUSE_NBHALO=OFF \
    -B build
cmake --build build -j
cmake --install build
(mkdir -p ./test &&
   cd ./test &&
   cp ../examples/template/param.in param.in &&
   cp ../examples/template/gptlnl gptlnl &&
   mpirun -n 4 ../build/install/bin/PowerLLEL 2>&1 | tee test.log)


###########################################
# Compile command for dependency software #
###########################################

#  For FFTW
# $HOME/fgn/software/fftw/3.3.8/configure \
#     --prefix=$HOME/fgn/software/fftw/3.3.8/build/install \
#     CC=icc \
#     CFLAGS="-O3 -xHost" \
#     MPICC=mpicc \
#     F77=ifort \
#     FFLAGS="-O3 -xHost" \
#     MAKEINFO=true \
#     --enable-avx512 \
#     --enable-openmp \
#     --enable-mpi

# For zlib
# cmake \
#     -DCMAKE_INSTALL_PREFIX=`pwd`/install \
#     -DCMAKE_C_COMPILER=icc \
#     -DCMAKE_C_FLAGS="-O3 -xHost" \
#     -LH \
#     ..

# For HDF5
# $HOME/fgn/software/hdf5/1.10.4/configure \
#     --with-zlib=$HOME/fgn/software/zlib/1.2.11/build/install \
#     --prefix=$HOME/fgn/software/hdf5/1.10.4/build/install \
#     --enable-fortran \
#     --enable-parallel \
#     CC=mpicc \
#     CFLAGS="-O3 -xHost" \
#     FC=mpifort \
#     FFLAGS="-O3 -xHost "

# For GPTL
# $HOME/fgn/software/gptl/8.1.1/configure \
#     --prefix=$HOME/fgn/software/gptl/8.1.1/build/install \
#     CC=mpicc \
#     CFLAGS="-O3 -xHost" \
#     MPICC=mpicc \
#     FC=mpifort \
#     FCLAGS="-O3 -xHost" \
#     MAKEINFO=true \
#     --enable-openmp \
#     --enable-pmpi