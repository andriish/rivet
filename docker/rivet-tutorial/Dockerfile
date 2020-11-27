ARG RIVET_VERSION
ARG PYTHIA_VERSION
FROM hepstore/rivet-pythia:${RIVET_VERSION}-${PYTHIA_VERSION}
LABEL maintainer="rivet-developers@cern.ch"

CMD /bin/bash

WORKDIR /work

ADD . /work
RUN rm /work/Dockerfile
