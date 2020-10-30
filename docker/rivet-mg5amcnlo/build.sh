#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.2
MG5_VERSION=2.7.3

BUILD="docker build ."
test "$TEST" && BUILD="echo $BUILD"

tag="hepstore/rivet-mg5amcnlo:${RIVET_VERSION}-${MG5_VERSION}"
echo "Building $tag"
$BUILD -f Dockerfile -t $tag

docker tag $tag hepstore/rivet-mg5amcnlo:$RIVET_VERSION
docker tag $tag hepstore/rivet-mg5amcnlo:latest

if [[ "$PUSH" = 1 ]]; then
    docker push $tag
    docker push hepstore/rivet-mg5amcnlo:latest
    docker push hepstore/rivet-mg5amcnlo:$RIVET_VERSION
fi
