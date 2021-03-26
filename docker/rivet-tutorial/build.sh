#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.3
PYTHIA_VERSION=8303

BUILD="docker build ."
BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILD="$BUILD --build-arg PYTHIA_VERSION=${PYTHIA_VERSION}"
test "$TEST" && BUILD="echo $BUILD"

tag="hepstore/rivet-tutorial:${RIVET_VERSION}"
echo "Building $tag"
$BUILD -f Dockerfile -t $tag

if [[ "$LATEST" = 1 ]]; then
    docker tag $tag hepstore/rivet-tutorial:latest
fi

if [[ "$PUSH" = 1 ]]; then
    docker push $tag
    if [[ "$LATEST" = 1 ]]; then
        sleep 1m
        docker push hepstore/rivet-tutorial:latest
    fi
fi
