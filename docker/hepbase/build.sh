#! /usr/bin/env bash

set -e

BUILD="docker build" # --squash"
test "$TEST" && BUILD="echo $BUILD"

for vhepmc in 2 3; do
    for tex in 0 1; do
    # for tex in 0; do

        MSG="Building hepbase image with HepMC=$vhepmc and TeX=$tex"

        BASEARGS="--build-arg RIVET_VERSION=3.1.3 --build-arg HEPMC_VERSION=$vhepmc --build-arg LATEX=$tex"
        GCCARGS="$BASEARGS --build-arg CXX_CMD=g++ --build-arg CC_CMD=gcc --build-arg FC_CMD=gfortran"
        CLANGARGS="$BASEARGS --build-arg CXX_CMD=clang++ --build-arg CC_CMD=clang --build-arg FC_CMD=gfortran"
        if [[ "$tex" = 1 ]]; then TEXSUFFIX="-latex"; else TEXSUFFIX=""; fi

        echo "@@ $MSG on Ubuntu with GCC compilers"
        tag=hepstore/hepbase-ubuntu-gcc-hepmc${vhepmc}-py3$TEXSUFFIX
        $BUILD . -f Dockerfile.ubuntu $GCCARGS -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        echo "@@ $MSG on Fedora with GCC compilers"
        tag=hepstore/hepbase-fedora-gcc-hepmc${vhepmc}-py3$TEXSUFFIX
        $BUILD . -f Dockerfile.fedora $GCCARGS -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        echo "@@ $MSG on CentOS with GCC compilers"
        tag=hepstore/hepbase-centos-gcc-hepmc${vhepmc}-py3$TEXSUFFIX
        $BUILD . -f Dockerfile.centos $GCCARGS -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        echo "@@ $MSG on Ubuntu with clang compilers"
        tag=hepstore/hepbase-ubuntu-clang-hepmc${vhepmc}-py3$TEXSUFFIX
        $BUILD . -f Dockerfile.ubuntu ${CLANGARGS/gfortran/flang} -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        # echo "@@ $MSG on Fedora with clang/LLVM compilers"
        # tag=hepstore/hepbase-fedora-clang-hepmc${vhepmc}-py3$TEXSUFFIX
        # $BUILD . -f Dockerfile.fedora ${CLANGARGS} -t $tag
        # test "$PUSH" = 1 && docker push $tag && sleep 1m
        # echo -e "\n\n\n"

        echo "@@ $MSG on Ubuntu with GCC compilers and Python 2"
        tag=hepstore/hepbase-ubuntu-gcc-hepmc${vhepmc}-py2$TEXSUFFIX
        $BUILD . -f Dockerfile.ubuntu-py2 $GCCARGS -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        # echo "@@ $MSG on Debian with GCC compilers and Python 2"
        # tag=hepstore/hepbase-debian-gcc-hepmc${vhepmc}-py2$TEXSUFFIX
        # $BUILD . -f Dockerfile.debpy2 $GCCARGS -t $tag
        # test "$PUSH" = 1 && docker push $tag && sleep 1m
        # echo -e "\n\n\n"

    done
done
