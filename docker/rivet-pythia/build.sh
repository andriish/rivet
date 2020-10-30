#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.2
PY8_VERSION=8303

BUILD="docker build ."
test "$TEST" && BUILD="echo $BUILD"

tag="hepstore/rivet-pythia:${RIVET_VERSION}-${PY8_VERSION}"
echo "Building $tag"
$BUILD -f Dockerfile -t $tag

docker tag $tag hepstore/rivet-pythia:$RIVET_VERSION
docker tag $tag hepstore/rivet-pythia:latest

if [[ "$PUSH" = 1 ]]; then
    docker push $tag
    sleep 1m
    docker push hepstore/rivet-pythia:latest
    sleep 1m
    docker push hepstore/rivet-pythia:$RIVET_VERSION
fi
