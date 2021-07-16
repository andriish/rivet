ARG RIVET_VERSION
FROM hepstore/rivet:${RIVET_VERSION}
LABEL maintainer="rivet-developers@cern.ch"
SHELL ["/bin/bash", "-c"]

ARG PYTHIA_VERSION

#COPY main93.cc Pythia8Rivet.h /
RUN export DEBIAN_FRONTEND=noninteractive \
    && mkdir /code && cd /code \
    && wget --no-verbose https://pythia.org/download/pythia83/pythia${PYTHIA_VERSION}.tgz -O- | tar xz

RUN cd /code \
    && cd pythia*/ && ./configure --{prefix,with-{hepmc2,lhapdf6,rivet}}=/usr/local \
    && make -j $(nproc --ignore=1) && make install \
    && cd examples \
    && make main93 && cp main93 /usr/local/bin/pythia8-main93 \
    && make main300 && cp main300 /usr/local/bin/pythia8-main300 \
    && rm -r /code

#&& mv /main93.cc . \
#&& mv /Pythia8Rivet.h include/Pythia8Plugins/ \

WORKDIR /work
