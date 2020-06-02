#! /usr/bin/env bash

set -e

for vhepmc in 2 3; do
    MSG="Building hepbase image with HepMC=$vhepmc"
    BASEARGS="--build-arg HEPMC_VERSION=$vhepmc"
    GCCARGS="$BASEARGS --build-arg CXX_CMD=g++ --build-arg CC_CMD=gcc --build-arg FC_CMD=gfortran"
    CLANGARGS="$BASEARGS --build-arg CXX_CMD=clang++ --build-arg CC_CMD=clang --build-arg FC_CMD=gfortran"

    echo "@@ $MSG on Ubuntu with GCC compilers"
    docker build . -f Dockerfile.ubuntu $GCCARGS -t hepstore/hepbase-ubuntu-gcc-hepmc${vhepmc}-py3

    echo "@@ $MSG on Fedora with GCC compilers"
    docker build . -f Dockerfile.fedora $GCCARGS -t hepstore/hepbase-fedora-gcc-hepmc${vhepmc}-py3

    echo "@@ $MSG on Ubuntu with clang compilers"
    docker build . -f Dockerfile.ubuntu ${CLANGARGS/gfortran/flang} -t hepstore/hepbase-ubuntu-clang-hepmc${vhepmc}-py3

    # echo "@@ $MSG on Fedora with clang/LLVM compilers"
    # docker build . -f Dockerfile.fedora ${CLANGARGS}                -t hepstore/hepbase-fedora-clang-hepmc${vhepmc}-py3

    echo "@@ $MSG on Debian with GCC compilers and Python 2"
    docker build . -f Dockerfile.debpy2 $GCCARGS -t hepstore/hepbase-debian-gcc-hepmc${vhepmc}-py2
done
