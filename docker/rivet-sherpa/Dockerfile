ARG RIVET_VERSION
FROM hepstore/rivet:${RIVET_VERSION}-ubuntu-gcc-hepmc2-py2
LABEL maintainer="rivet-developers@cern.ch"

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y rsync pypy \
    && cd /usr/local \
    && git clone https://gitlab.com/openloops/OpenLoops.git \
    && cd OpenLoops && ./scons \
    && ./openloops libinstall ppll pptt

ARG SHERPA_VERSION
ENV PATH="${PATH}:/usr/lib64/openmpi/bin"
RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y libtool texinfo libsqlite3-dev swig openmpi-bin libopenmpi-dev \
    && mkdir /code && cd /code \
    && git clone -b v${SHERPA_VERSION} https://gitlab.com/sherpa-team/sherpa.git \
    && cd sherpa && autoreconf -fi \
    && ./configure --prefix=/usr/local \
      --enable-lhapdf=/usr/local --enable-fastjet=/usr/local \
      --enable-hepmc2=/usr/local --enable-openloops=/usr/local/OpenLoops \
      --enable-rivet=/usr/local --enable-analysis \
      --enable-pyext --enable-gzip \
      --enable-mpi CC=mpicc CXX=mpic++ \
    && make -j $(nproc --ignore=1) CXXFLAGS="-O2" AM_CXXFLAGS="-O2" \
    && make install \
    && mkdir /usr/local/Sherpa \
    && mv Examples /usr/local/Sherpa \
    && cd / && rm -r /code


# RUN apt-get update; \
# apt-get -q -y --no-install-recommends install pypy; \
# apt autoremove; apt-get autoclean; \
# rm -rf /var/lib/apt/lists /var/cache/apt; \
# rm -rf $(find / -name doc 2>&1 | grep -v Permission)

# RUN echo "export PATH=/usr/local/bin:\$PATH\n\
# export LD_LIBRARY_PATH=/usr/local/lib:\$LD_LIBRARY_PATH\n\
# . /usr/local/rivetenv.sh\n\
# if [ "\$PS1" ]; then PS1='\\h:\\w\\$ '; fi" > /etc/bash.bashrc

WORKDIR /work
