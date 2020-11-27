#! /usr/bin/env bash

set -e

VERSION=3.1.3
BUILD="docker build . -f Dockerfile"
BUILD="$BUILD --build-arg RIVET_VERSION=$VERSION" # --squash"
test "$TEST" && BUILD="echo $BUILD"

MSG="Building Rivet $VERSION image with architecture ="
for vhepmc in 3 2; do
    for arch in ubuntu-gcc-hepmc$vhepmc-py3 ubuntu-gcc-hepmc$vhepmc-py3 ubuntu-gcc-hepmc$vhepmc-py2 ubuntu-clang-hepmc$vhepmc-py3; do

        echo "@@ $MSG $arch"
        tag="hepstore/rivet:$VERSION-$arch"
        $BUILD --build-arg ARCH=$arch -t $tag
        if [[ "$PUSH" = 1 ]]; then
            docker push $tag
            sleep 1m
        fi
        echo -e "\n\n\n"

    done
done

## Convenience tags
docker tag hepstore/rivet:$VERSION{-ubuntu-gcc-hepmc2-py3,}
docker tag hepstore/rivet:$VERSION{-ubuntu-gcc-hepmc3-py3,-hepmc3}
if [[ "$LATEST" = 1 ]]; then
    docker tag hepstore/rivet:{$VERSION,latest}
fi
if [[ "$PUSH" = 1 ]]; then
    docker push hepstore/rivet:$VERSION
    sleep 1m
    docker push hepstore/rivet:$VERSION-hepmc3
    if [[ "$LATEST" = 1 ]]; then
        sleep 1m
        docker push hepstore/rivet:latest
    fi
fi
