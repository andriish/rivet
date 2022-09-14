#! /usr/bin/env bash

set -e

BUILD="docker build . -f Dockerfile"
function xdocker { echo "docker $@"; docker "$@"; }

test "$FORCE" && BUILD="$BUILD --no-cache"

YODA_BRANCH=yoda-1.9.6

## Last branch name -> latest
for RIVET_BRANCH in release-3-1-x rivet-3.1.6; do
    RIVET_VERSION=${RIVET_BRANCH#rivet-}

    BUILD="$BUILD --build-arg YODA_BRANCH=$YODA_BRANCH" # --squash"
    BUILD="$BUILD --build-arg RIVET_BRANCH=$RIVET_BRANCH" # --squash"
    test "$TEST" = 1 && BUILD="echo $BUILD"

    test "$INTEL" = 1 && intel="ubuntu-intel-hepmc$vhepmc-py3"

    MSG="Building Rivet $RIVET_VERSION image with architecture ="
    for vhepmc in 3; do   # 2
        for arch in $intel ubuntu-gcc-hepmc$vhepmc-py3 ubuntu-gcc-hepmc$vhepmc-py3 ubuntu-clang-hepmc$vhepmc-py3 ubuntu-gcc-hepmc$vhepmc-py2; do

            echo "@@ $MSG $arch"
            tag="hepstore/rivet:$RIVET_VERSION-$arch"
            $BUILD --build-arg ARCH=$arch -t $tag
            if [[ "$PUSH" = 1 ]]; then
                xdocker push $tag
                test "$NOSLEEP" = 1 || sleep 1m
            fi
            echo -e "\n\n\n"

        done
    done
done

## Convenience tags
xdocker tag hepstore/rivet:$RIVET_VERSION{-ubuntu-gcc-hepmc3-py3,-hepmc3}
xdocker tag hepstore/rivet:$RIVET_VERSION{-hepmc3,}
if [[ "$LATEST" = 1 ]]; then
    xdocker tag hepstore/rivet:{$RIVET_VERSION,latest}
fi
if [[ "$PUSH" = 1 ]]; then
    xdocker push hepstore/rivet:$RIVET_VERSION
    test "$NOSLEEP" = 1 || sleep 1m
    xdocker push hepstore/rivet:$RIVET_VERSION-hepmc3
    if [[ "$LATEST" = 1 ]]; then
        test "$NOSLEEP" = 1 || sleep 1m
        xdocker push hepstore/rivet:latest
    fi
fi
