#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.2

BUILD="docker build ."
test "$TEST" && BUILD="echo $BUILD"

tag="hepstore/rivet-tutorial:${RIVET_VERSION}"
echo "Building $tag"
$BUILD -f Dockerfile -t $tag

docker tag $tag hepstore/rivet-tutorial:latest

if [[ "$PUSH" = 1 ]]; then
    docker push $tag
    sleep 1m
    docker push hepstore/rivet-tutorial:latest
fi
