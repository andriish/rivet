ARG RIVET_VERSION
FROM hepstore/rivet:${RIVET_VERSION}
LABEL maintainer="rivet-developers@cern.ch"
SHELL ["/bin/bash", "-c"]

ARG PYTHIA_VERSION

COPY main89R.cc Makefile.main89R /code/
RUN export DEBIAN_FRONTEND=noninteractive \
    && cd /code \
    && wget --no-verbose https://pythia.org/download/pythia83/pythia${PYTHIA_VERSION}.tgz -O- | tar xz \
    && cd pythia*/ && ./configure --{prefix,with-{hepmc3,lhapdf6,rivet}}=/usr/local \
    && cd examples \
    && cp /code/main89R.cc ./ \
    && cat /code/Makefile.main89R >> Makefile \
    && cd .. \
    && make -j $(nproc --ignore=1) && make install \
    && cd examples \
    && make main89R && cp main89R /usr/local/bin/pythia8-main89R \
    && make main93 && cp main93 /usr/local/bin/pythia8-main93 \
    && make main300 && cp main300 /usr/local/bin/pythia8-main300 \
    && rm -r /code

WORKDIR /work
