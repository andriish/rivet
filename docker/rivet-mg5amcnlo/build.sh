#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.4
MG5_VERSION=2.7.3   #< also edit in Dockerfile, for now

BUILD="docker build ."
BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION}"
test "$TEST" && BUILD="echo $BUILD"

tag="hepstore/rivet-mg5amcnlo:${RIVET_VERSION}-${MG5_VERSION}"
echo "Building $tag"
$BUILD -f Dockerfile -t $tag

docker tag $tag hepstore/rivet-mg5amcnlo:$RIVET_VERSION
if [[ "$LATEST" = 1 ]]; then
    docker tag $tag hepstore/rivet-mg5amcnlo:latest
fi

if [[ "$PUSH" = 1 ]]; then
    docker push $tag
    sleep 1m
    docker push hepstore/rivet-mg5amcnlo:$RIVET_VERSION
    if [[ "$LATEST" = 1 ]]; then
        sleep 1m
        docker push hepstore/rivet-mg5amcnlo:latest
    fi
fi
