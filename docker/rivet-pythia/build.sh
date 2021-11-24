#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.5
PYTHIA_VERSION=8306

BUILD="docker build ."

test "$FORCE" && BUILD="$BUILD --force-rm"

BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILD="$BUILD --build-arg PYTHIA_VERSION=${PYTHIA_VERSION}"
test "$TEST" && BUILD="echo $BUILD"

tag="hepstore/rivet-pythia:${RIVET_VERSION}-${PYTHIA_VERSION}"
echo "Building $tag"
$BUILD -f Dockerfile -t $tag

docker tag $tag hepstore/rivet-pythia:$RIVET_VERSION
if [[ "$LATEST" = 1 ]]; then
    docker tag $tag hepstore/rivet-pythia:latest
fi

if [[ "$PUSH" = 1 ]]; then
    docker push $tag
    sleep 1m
    docker push hepstore/rivet-pythia:$RIVET_VERSION
    if [[ "$LATEST" = 1 ]]; then
        sleep 1m
        docker push hepstore/rivet-pythia:latest
    fi
fi
