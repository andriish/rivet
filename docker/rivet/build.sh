#! /usr/bin/env bash

set -e

VERSION=3.1.2
BUILD="docker build . -f Dockerfile --build-arg RIVET_VERSION=$VERSION" # --squash"
test "$TEST" && BUILD="echo $BUILD"

MSG="Building Rivet $VERSION image with architecture ="
for vhepmc in 3 2; do
    for arch in ubuntu-gcc-hepmc$vhepmc-py3 ubuntu-gcc-hepmc$vhepmc-py3 ubuntu-gcc-hepmc$vhepmc-py2 ubuntu-clang-hepmc$vhepmc-py3; do

        echo "@@ $MSG $arch"
        tag="hepstore/rivet:$VERSION-$arch"
        $BUILD --build-arg ARCH=$arch -t $tag
        test "$PUSH" = 1 && docker push $tag && sleep 1m
        echo -e "\n\n\n"

    done
done

## Convenience tags
docker tag hepstore/rivet:$VERSION{-ubuntu-gcc-hepmc2-py3,}
docker tag hepstore/rivet:$VERSION{-ubuntu-gcc-hepmc3-py3,-hepmc3}
docker tag hepstore/rivet:{$VERSION,latest}
if [[ "$PUSH" = 1 ]]; then
    docker push hepstore/rivet:$VERSION
    sleep 30s
    docker push hepstore/rivet:$VERSION-hepmc3
    sleep 30s
    docker push hepstore/rivet:latest
fi
