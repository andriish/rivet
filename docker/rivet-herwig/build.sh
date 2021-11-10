#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.5
YODA_VERSION=1.9.2
THEPEG_VERSION=2.2.2
HERWIG_VERSION=7.2.2

BUILD="docker build ."
BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILD="$BUILD --build-arg YODA_VERSION=${YODA_VERSION}"
BUILD="$BUILD --build-arg THEPEG_VERSION=${THEPEG_VERSION}"
BUILD="$BUILD --build-arg HERWIG_VERSION=${HERWIG_VERSION}"
test "$TEST" && BUILD="echo $BUILD"

tag3a="hepstore/rivet-herwig:${RIVET_VERSION}-${HERWIG_VERSION}"
tag3b="$tag3a-py3"
echo "Building $tag3a"
$BUILD -f Dockerfile.ubuntu -t $tag3a
docker tag $tag3a $tag3b

echo -e "\n\n"

# tag2="hepstore/rivet-herwig:${RIVET_VERSION}-${HERWIG_VERSION}-py2"
# echo "Building $tag2"
# $BUILD -f Dockerfile.ubuntu-py2 -t $tag2

docker tag $tag3b hepstore/rivet-herwig:$RIVET_VERSION
if [[ "$LATEST" = 1 ]]; then
    docker tag $tag3b hepstore/rivet-herwig:latest
fi

if [[ "$PUSH" = 1 ]]; then
    docker push $tag3a
    sleep 30s
    docker push $tag3b
    # sleep 1m
    # docker push $tag2
    sleep 1m
    docker push hepstore/rivet-herwig:$RIVET_VERSION
    if [[ "$LATEST" = 1 ]]; then
        sleep 1m
        docker push hepstore/rivet-herwig:latest
    fi
fi
