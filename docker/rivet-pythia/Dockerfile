FROM hepstore/rivet:3.1.2
LABEL maintainer="rivet-developers@cern.ch"
SHELL ["/bin/bash", "-c"]

ARG PY8_VERSION=8303

# TODO: remove with next hepbase/rivet image version
RUN export DEBIAN_FRONTEND=noninteractive \
    && mkdir /code && cd /code \
    && wget --no-verbose https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.3.0.tar.gz -O- | tar xz \
    && cd LHAPDF-*/ && ./configure --prefix=/usr/local \
    && make -j $(nproc --ignore=1) && make install \
    && rm -r /code

# RUN export DEBIAN_FRONTEND=noninteractive \
#     && mkdir /code && cd /code \
#     && apt-get -y install rsync && apt-get autoclean \
#     && wget --no-verbose http://home.thep.lu.se/~torbjorn/pythia8/pythia${PY8_VERSION}.tgz -O- | tar xz \
#     && cd pythia*/ && ./configure --enable-shared --{prefix,with-{hepmc2,lhapdf6,rivet}}=/usr/local \
#     && make -j $(nproc --ignore=1) && make install \
#     && cd examples \
#     && make main93 && cp main93 /usr/local/bin/pythia8-main93 \
#     && make main300 && cp main300 /usr/local/bin/pythia8-main300 \
#     && rm -r /code

COPY main93.cc Pythia8Rivet.h /
RUN export DEBIAN_FRONTEND=noninteractive \
    && mkdir /code && cd /code \
    && apt-get -y install rsync && apt-get autoclean \
    && wget --no-verbose http://home.thep.lu.se/~torbjorn/pythia8/pythia${PY8_VERSION}.tgz -O- | tar xz \
    && cd pythia*/ && ./configure --enable-shared --{prefix,with-{hepmc2,lhapdf6,rivet}}=/usr/local \
    && mv /Pythia8Rivet.h include/Pythia8Plugins/ \
    && make -j $(nproc --ignore=1) && make install \
    && cd examples \
    && mv /main93.cc . \
    && make main93 && cp main93 /usr/local/bin/pythia8-main93 \
    && make main300 && cp main300 /usr/local/bin/pythia8-main300 \
    && rm -r /code

# COPY main93.cc /code
# RUN export DEBIAN_FRONTEND=noninteractive \
#     && cd /code \
#     && cp main93.cc pythia*/examples/ \
#     && cd pythia*/examples \
#     && make main93 && cp main93 /usr/local/bin/pythia8-main93

# RUN export DEBIAN_FRONTEND=noninteractive \
#     && mkdir /code && cd /code \
#     && wget --no-verbose https://agile.hepforge.org/downloads/?f=Sacrifice-1.1.2.tar.gz -O- | tar xz \
#     && cd Sacrifice-*/ \
#     && ./configure --prefix=/usr/local --with-pythia=/usr/local --with-hepmc=/usr/local --with-lhapdf=/usr/local \
#     && make -j $(nproc --ignore=1) && make install \
#     && cd ../.. && rm -r /code

WORKDIR /work
