#! /usr/bin/env bash

# for vpython in 2 3; do
#     for vhepmc in 2.06.10 3.2.1; do
#         vhepmc_major=$(echo $vhepmc | cut -d. -f1)
#         for cxx in gcc-g++ clang; do
#             compiler=$(echo $cxx | cut -d- -f1)
#             docker build . -f Dockerfile.fedora --build-arg PYTHON_VERSION=$vpython --build-arg HEPMC_VERSION=$vhepmc --build-arg CXX_COMPILER=$cxx -t hepstore/fedora-${compiler}-py${vpython}-hepmc${vhepmc_major}
#         done
#     done
# done

VPYTHONS="2" # 3"
VHEPMCS="2" # 3"

for vpython in $VPYTHONS; do
    for vhepmc in $VHEPMCS; do
        MSG="Building hepbase image with Python=$vpython, HepMC=$vhepmc"
        BASEARGS="--build-arg PYTHON_VERSION=$vpython --build-arg HEPMC_VERSION=$vhepmc"

        GCCARGS="$BASEARGS --build-arg CXX_CMD=g++ --build-arg CC_CMD=gcc --build-arg FC=gfortran"
        echo "@@ $MSG on Ubuntu with GCC compilers"
        docker build . -f Dockerfile.ubuntu $GCCARGS -t hepstore/ubuntu-gcc-py${vpython}-hepmc${vhepmc}
        echo "@@ $MSG on Fedora with GCC compilers"
        docker build . -f Dockerfile.fedora $GCCARGS -t hepstore/fedora-gcc-py${vpython}-hepmc${vhepmc}

        CLANGARGS="$BASEARGS --build-arg CXX_CMD=clang++ --build-arg CC_CMD=clang --build-arg FC=gfortran"
        echo "@@ $MSG on Ubuntu with clang compilers"
        docker build . -f Dockerfile.ubuntu ${CLANGARGS/gfortran/flang} -t hepstore/ubuntu-clang-py${vpython}-hepmc${vhepmc}
        echo "@@ $MSG on Fedora with clang/LLVM compilers"
        docker build . -f Dockerfile.fedora ${CLANGARGS}                -t hepstore/fedora-clang-py${vpython}-hepmc${vhepmc}
    done
done
