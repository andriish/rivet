FROM hepstore/rivet:3.1.2-ubuntu-gcc-hepmc2-py2

LABEL maintainer="rivet-developers@cern.ch"

# TODO: remove with next hepbase/rivet image version
RUN export DEBIAN_FRONTEND=noninteractive \
    && mkdir /code && cd /code \
    && wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.3.0.tar.gz -O- | tar xz \
    && cd LHAPDF-*/ && ./configure --prefix=/usr/local \
    && make -j $(nproc --ignore=1) && make install \
    && cd ../.. && rm -r /code

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y rsync pypy \
    && cd /usr/local \
    && git clone https://gitlab.com/openloops/OpenLoops.git \
    && cd OpenLoops && ./scons \
    && ./openloops libinstall ppll pptt

ENV PATH="${PATH}:/usr/lib64/openmpi/bin"
RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y libtool texinfo libsqlite3-dev swig openmpi-bin libopenmpi-dev \
    && mkdir /code && cd /code \
    && git clone -b rel-2-2-10 https://gitlab.com/sherpa-team/sherpa.git \
    && cd sherpa && autoreconf -fi \
    && ./configure --prefix=/usr/local \
      --enable-lhapdf=/usr/local --enable-fastjet=/usr/local \
      --enable-hepmc2=/usr/local --enable-openloops=/usr/local/OpenLoops \
      --enable-pyext --enable-gzip --enable-analysis \
      --enable-mpi CC=mpicc CXX=mpic++ \
    && make -j $(nproc --ignore=1) CXXFLAGS="-O2 -std=c++11" AM_CXXFLAGS="-O2 -std=c++11" \
    && make install \
    && mkdir /usr/local/Sherpa \
    && mv Examples /usr/local/Sherpa


# && cd ../.. && rm -r /code

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
