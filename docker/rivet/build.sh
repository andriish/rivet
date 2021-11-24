#! /usr/bin/env bash

set -e

BUILD="docker build . -f Dockerfile"

test "$FORCE" && BUILD="$BUILD --force-rm"

## Last branch name -> latest
for RIVET_BRANCH in release-3-1-x rivet-3.1.5; do
    RIVET_VERSION=${RIVET_BRANCH#rivet-}

    BUILD="$BUILD --build-arg RIVET_BRANCH=$RIVET_BRANCH" # --squash"
    test "$TEST" && BUILD="echo $BUILD"

    MSG="Building Rivet $RIVET_VERSION image with architecture ="
    for vhepmc in 3; do   # 2
        for arch in ubuntu-intel-hepmc$vhepmc-py3 ubuntu-gcc-hepmc$vhepmc-py3 ubuntu-gcc-hepmc$vhepmc-py3 ubuntu-clang-hepmc$vhepmc-py3 ubuntu-gcc-hepmc$vhepmc-py2; do

            echo "@@ $MSG $arch"
            tag="hepstore/rivet:$RIVET_VERSION-$arch"
            $BUILD --build-arg ARCH=$arch -t $tag
            if [[ "$PUSH" = 1 ]]; then
                docker push $tag
                sleep 1m
            fi
            echo -e "\n\n\n"

        done
    done
done

## Convenience tags
docker tag hepstore/rivet:$RIVET_VERSION{-ubuntu-gcc-hepmc3-py3,-hepmc3}
docker tag hepstore/rivet:$RIVET_VERSION{-hepmc3,}
if [[ "$LATEST" = 1 ]]; then
    docker tag hepstore/rivet:{$RIVET_VERSION,latest}
fi
if [[ "$PUSH" = 1 ]]; then
    docker push hepstore/rivet:$RIVET_VERSION
    sleep 1m
    docker push hepstore/rivet:$RIVET_VERSION-hepmc3
    if [[ "$LATEST" = 1 ]]; then
        sleep 1m
        docker push hepstore/rivet:latest
    fi
fi
