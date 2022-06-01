#! /usr/bin/env bash

set -e

BUILD="docker build" # --squash"

test "$FORCE" && BUILD="$BUILD --no-cache"

test "$TEST" && BUILD="echo $BUILD"

# FEDORA_IMAGE=fedora:32
# UBUNTU_IMAGE=ubuntu:22.04

RIVETBS_VERSION=3.1.6
LHAPDF_VERSION=6.5.1

for vhepmc in 3.2.5; do   # 2.06.11
    for tex in 0 1; do

        MSG="Building hepbase image with HepMC=$vhepmc and TeX=$tex"

        BASEARGS="--build-arg RIVETBS_VERSION=$RIVETBS_VERSION"
        BASEARGS="$BASEARGS --build-arg HEPMC_VERSION=$vhepmc"
        BASEARGS="$BASEARGS --build-arg LHAPDF_VERSION=$LHAPDF_VERSION"
        BASEARGS="$BASEARGS --build-arg LATEX=$tex"

        GCCARGS="$BASEARGS --build-arg BUILD_TOOLS=GCC"
        CLANGARGS="$BASEARGS --build-arg BUILD_TOOLS=LLVM"
        INTELARGS="$BASEARGS --build-arg BUILD_TOOLS=Intel"
        if [[ "$tex" = 1 ]]; then TEXSUFFIX="-latex"; else TEXSUFFIX=""; fi

        echo "@@ $MSG on Ubuntu with GCC compilers"
        tag=hepstore/hepbase-ubuntu-gcc-hepmc${vhepmc:0:1}-py3$TEXSUFFIX
        $BUILD . -f Dockerfile.ubuntu $GCCARGS -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        echo "@@ $MSG on Ubuntu LTS 20.04 with GCC compilers"
        tag=hepstore/hepbase-ubuntu20-gcc-hepmc${vhepmc:0:1}-py3$TEXSUFFIX
        $BUILD . -f Dockerfile.ubuntu $GCCARGS -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        echo "@@ $MSG on Ubuntu with clang compilers"
        tag=hepstore/hepbase-ubuntu-clang-hepmc${vhepmc:0:1}-py3$TEXSUFFIX
        $BUILD . -f Dockerfile.ubuntu ${CLANGARGS/gfortran/flang} -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        if [[ "$INTEL" = 1 ]]; then
            echo "@@ $MSG on Ubuntu with Intel compilers"
            tag=hepstore/hepbase-ubuntu-intel-hepmc${vhepmc:0:1}-py3$TEXSUFFIX
            $BUILD . -f Dockerfile.ubuntu $INTELARGS -t $tag
            test "$PUSH" = 1 && docker push $tag && sleep 1m
            echo -e "\n\n\n"
        fi

        echo "@@ $MSG on Fedora with GCC compilers"
        tag=hepstore/hepbase-fedora-gcc-hepmc${vhepmc:0:1}-py3$TEXSUFFIX
        $BUILD . -f Dockerfile.fedora $GCCARGS -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

        # echo "@@ $MSG on Ubuntu rolling with GCC compilers"
        # tag=hepstore/hepbase-ubuntux-gcc-hepmc${vhepmc:0:1}-py3$TEXSUFFIX
        # $BUILD . -f Dockerfile.ubuntu $GCCARGS --build-arg UBUNTU_IMAGE=ubuntu:rolling -t $tag
        # test "$PUSH" = 1 && docker push $tag && sleep 1m
        # echo -e "\n\n\n"

        # echo "@@ $MSG on Fedora rawhide with GCC compilers"
        # tag=hepstore/hepbase-fedorax-gcc-hepmc${vhepmc:0:1}-py3$TEXSUFFIX
        # $BUILD . -f Dockerfile.fedora $GCCARGS --build-arg FEDORA_IMAGE=fedora:rawhide -t $tag
        # test "$PUSH" = 1 && docker push $tag && sleep 1m
        # echo -e "\n\n\n"

        echo "@@ $MSG on Ubuntu with GCC compilers and Python 2"
        tag=hepstore/hepbase-ubuntu-gcc-hepmc${vhepmc:0:1}-py2$TEXSUFFIX
        $BUILD . -f Dockerfile.ubuntu-py2 $GCCARGS -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

    done
done
