#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.6
YODA_VERSION=1.9.5
THEPEG_VERSION=2.2.3
HERWIG_VERSION=7.2.3
LHAPDF_VERSION=6.5.2

BUILD="docker build ."

test "$FORCE" && BUILD="$BUILD --no-cache"

BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILD="$BUILD --build-arg YODA_VERSION=${YODA_VERSION}"
BUILD="$BUILD --build-arg THEPEG_VERSION=${THEPEG_VERSION}"
BUILD="$BUILD --build-arg HERWIG_VERSION=${HERWIG_VERSION}"
BUILD="$BUILD --build-arg LHAPDF_VERSION=${LHAPDF_VERSION}"
test "$TEST" && BUILD="echo $BUILD"

tag3a="hepstore/rivet-herwig:${RIVET_VERSION}-${HERWIG_VERSION}"
tag3b="$tag3a-py3"
echo "Building $tag3a"
$BUILD -f Dockerfile -t $tag3a
docker tag $tag3a $tag3b

echo -e "\n\n"

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
