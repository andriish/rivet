#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.6
SHERPA_VERSION=2.2.12

BUILD="docker build ."

test "$FORCE" && BUILD="$BUILD --no-cache"

BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILD="$BUILD --build-arg SHERPA_VERSION=${SHERPA_VERSION}"
test "$TEST" && BUILD="echo $BUILD"

tag="hepstore/rivet-sherpa:${RIVET_VERSION}-${SHERPA_VERSION}"
echo "Building $tag"
$BUILD -f Dockerfile -t $tag

docker tag $tag hepstore/rivet-sherpa:$RIVET_VERSION
if [[ "$LATEST" = 1 ]]; then
    docker tag $tag hepstore/rivet-sherpa:latest
fi

if [[ "$PUSH" = 1 ]]; then
    docker push $tag
    sleep 1m
    docker push hepstore/rivet-sherpa:$RIVET_VERSION
    if [[ "$LATEST" = 1 ]]; then
        sleep 1m
        docker push hepstore/rivet-sherpa:latest
    fi
fi
