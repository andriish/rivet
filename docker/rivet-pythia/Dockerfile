ARG RIVET_VERSION
FROM hepstore/rivet:${RIVET_VERSION}
LABEL maintainer="rivet-developers@cern.ch"
SHELL ["/bin/bash", "-c"]

ARG PYTHIA_VERSION

COPY main93.cc Pythia8Rivet.h /
RUN export DEBIAN_FRONTEND=noninteractive \
    && mkdir /code && cd /code \
    && wget --no-verbose http://home.thep.lu.se/~torbjorn/pythia8/pythia${PYTHIA_VERSION}.tgz -O- | tar xz \
    && cd pythia*/ && ./configure --enable-shared --{prefix,with-{hepmc2,lhapdf6,rivet}}=/usr/local \
    && mv /Pythia8Rivet.h include/Pythia8Plugins/ \
    && make -j $(nproc --ignore=1) && make install \
    && cd examples \
    && mv /main93.cc . \
    && make main93 && cp main93 /usr/local/bin/pythia8-main93 \
    && make main300 && cp main300 /usr/local/bin/pythia8-main300 \
    && rm -r /code

WORKDIR /work
