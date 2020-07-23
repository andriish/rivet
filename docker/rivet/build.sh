#! /usr/bin/env bash

set -e

VERSION=3.1.2
BUILD="docker build . -f Dockerfile --build-arg RIVET_VERSION=$VERSION" # --squash"
test "$TEST" && BUILD="echo $BUILD"

MSG="Building Rivet $VERSION image with architecture ="
for vhepmc in 2 3; do
    for arch in ubuntu-gcc-hepmc$vhepmc-py3 ubuntu-gcc-hepmc$vhepmc-py3 ubuntu-gcc-hepmc$vhepmc-py2 ubuntu-clang-hepmc$vhepmc-py3; do

        echo "@@ $MSG $arch"
        tag="rivet:$VERSION-$arch"
        $BUILD --build-arg BASE=$arch -t $tag
        test "$PUSH" = 1 && docker push $tag
        echo -e "\n\n\n"

    done
done
