#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.2
YODA_VERSION=1.8.3
THEPEG_VERSION=2.2.1
HERWIG_VERSION=7.2.1

BUILD="docker build ."
BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION} --build-arg YODA_VERSION=${YODA_VERSION}"
BUILD="$BUILD --build-arg THEPEG_VERSION=${THEPEG_VERSION} --build-arg HERWIG_VERSION=${HERWIG_VERSION}"
test "$TEST" && BUILD="echo $BUILD"

tag3="hepstore/rivet-herwig:${RIVET_VERSION}-${HERWIG_VERSION}"
echo "Building $tag3"
$BUILD -f Dockerfile.ubuntu -t $tag3

echo -e "\n\n"

tag2="hepstore/rivet-herwig:${RIVET_VERSION}-${HERWIG_VERSION}-py2"
echo "Building $tag2"
$BUILD -f Dockerfile.ubuntu-py2 -t $tag2

docker tag $tag2 hepstore/rivet-herwig:$RIVET_VERSION
docker tag $tag2 hepstore/rivet-herwig:latest

if [[ "$PUSH" = 1 ]]; then
    docker push $tag3
    sleep 1m
    docker push $tag2
    sleep 1m
    docker push hepstore/rivet-herwig:$RIVET_VERSION
    sleep 1m
    docker push hepstore/rivet-herwig:latest
fi
