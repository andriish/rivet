#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.6
MG5_VERSION=3.2.0
MG5_URL=https://launchpad.net/mg5amcnlo/3.0/3.2.x/+download/MG5_aMC_v3.2.0.tar.gz

BUILD="docker build ."

test "$FORCE" && BUILD="$BUILD --no-cache"

BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION} --build-arg MG5_URL=${MG5_URL}"
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
