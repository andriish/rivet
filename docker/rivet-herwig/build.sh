#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.2
YODA_VERSION=1.8.3
THEPEG_VERSION=2.2.1
HERWIG_VERSION=7.2.1

BUILD="docker build . -f Dockerfile"
BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION} --build-arg YODA_VERSION=${YODA_VERSION}"
BUILD="$BUILD --build-arg THEPEG_VERSION=${THEPEG_VERSION} --build-arg HERWIG_VERSION=${HERWIG_VERSION}"
test "$TEST" && BUILD="echo $BUILD"

tag="hepstore/rivet-herwig:${RIVET_VERSION}-${HERWIG_VERSION}"
$BUILD -t $tag

docker tag $tag hepstore/rivet-herwig:$RIVET_VERSION
docker tag $tag hepstore/rivet-herwig:latest

if [[ "$PUSH" = 1 ]]; then
    docker push $tag
    docker push hepstore/rivet-herwig:$RIVET_VERSION
    docker push hepstore/rivet-herwig:latest
fi
