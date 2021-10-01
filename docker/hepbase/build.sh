#! /usr/bin/env bash

set -e

BUILD="docker build" # --squash"
test "$TEST" && BUILD="echo $BUILD"

for vhepmc in 3.2.4  2.06.11; do
    for tex in 0 1; do

        MSG="Building hepbase image with HepMC=$vhepmc and TeX=$tex"

        BASEARGS="--build-arg RIVETBS_VERSION=3.1.4 --build-arg HEPMC_VERSION=$vhepmc --build-arg LATEX=$tex"
        GCCARGS="$BASEARGS --build-arg BUILD_TOOLS=GCC"
        CLANGARGS="$BASEARGS --build-arg BUILD_TOOLS=LLVM"
        INTELARGS="$BASEARGS --build-arg BUILD_TOOLS=Intel"
        if [[ "$tex" = 1 ]]; then TEXSUFFIX="-latex"; else TEXSUFFIX=""; fi

        echo "@@ $MSG on Ubuntu with GCC compilers"
        tag=hepstore/hepbase-ubuntu-gcc-hepmc${vhepmc:0:1}-py3$TEXSUFFIX
        $BUILD . -f Dockerfile.ubuntu $GCCARGS -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        echo "@@ $MSG on Ubuntu with clang compilers"
        tag=hepstore/hepbase-ubuntu-clang-hepmc${vhepmc:0:1}-py3$TEXSUFFIX
        $BUILD . -f Dockerfile.ubuntu ${CLANGARGS/gfortran/flang} -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        echo "@@ $MSG on Ubuntu with Intel compilers"
        tag=hepstore/hepbase-ubuntu-intel-hepmc${vhepmc:0:1}-py3$TEXSUFFIX
        $BUILD . -f Dockerfile.ubuntu $INTELARGS -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        echo "@@ $MSG on Fedora with GCC compilers"
        tag=hepstore/hepbase-fedora-gcc-hepmc${vhepmc:0:1}-py3$TEXSUFFIX
        $BUILD . -f Dockerfile.fedora $GCCARGS -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        # echo "@@ $MSG on Ubuntu with GCC compilers and Python 2"
        # tag=hepstore/hepbase-ubuntu-gcc-hepmc${vhepmc:0:1}-py2$TEXSUFFIX
        # $BUILD . -f Dockerfile.ubuntu-py2 $GCCARGS -t $tag
        # test "$PUSH" = 1 && docker push $tag && sleep 1m
        # echo -e "\n\n\n"

    done
done
