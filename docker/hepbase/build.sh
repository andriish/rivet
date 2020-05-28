#! /usr/bin/env bash

set -e

VPYTHONS="2" # 3"
VHEPMCS="2" # 3"

for vpython in $VPYTHONS; do
    for vhepmc in $VHEPMCS; do
        MSG="Building hepbase image with Python=$vpython, HepMC=$vhepmc"
        BASEARGS="--build-arg PYTHON_VERSION=$vpython --build-arg HEPMC_VERSION=$vhepmc"

        # GCCARGS="$BASEARGS --build-arg CXX_CMD=g++ --build-arg CC_CMD=gcc --build-arg FC=gfortran"
        # echo "@@ $MSG on Ubuntu with GCC compilers"
        # docker build . -f Dockerfile.ubuntu $GCCARGS -t hepstore/ubuntu-gcc-py${vpython}-hepmc${vhepmc}
        # echo "@@ $MSG on Fedora with GCC compilers"
        # docker build . -f Dockerfile.fedora $GCCARGS -t hepstore/fedora-gcc-py${vpython}-hepmc${vhepmc}

        CLANGARGS="$BASEARGS --build-arg CXX_CMD=clang++ --build-arg CC_CMD=clang --build-arg FC=gfortran"
        # echo "@@ $MSG on Ubuntu with clang compilers"
        # docker build . -f Dockerfile.ubuntu ${CLANGARGS/gfortran/flang} -t hepstore/ubuntu-clang-py${vpython}-hepmc${vhepmc}
        echo "@@ $MSG on Fedora with clang/LLVM compilers"
        docker build . -f Dockerfile.fedora ${CLANGARGS}                -t hepstore/fedora-clang-py${vpython}-hepmc${vhepmc}
    done
done
