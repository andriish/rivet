FROM hepstore/rivet:3.1.1
LABEL maintainer="rivet-developers@cern.ch"

RUN mkdir /code && cd /code \
    && wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.2.3.tar.gz -O- | tar xz \
    && cd LHAPDF-*/ && ./configure --prefix=/usr/local \
    && make -j5 && make install \
    && cd ../.. && rm -r /code

RUN dnf install -y rsync git pypy \
    && cd /usr/local \
    && git clone https://gitlab.com/openloops/OpenLoops.git \
    && cd OpenLoops && ./scons \
    && ./openloops libinstall ppll pptt

RUN dnf install -y libtool texinfo sqlite-devel swig pypy openmpi openmpi-devel
ENV PATH="${PATH}:/usr/lib64/openmpi/bin"

RUN mkdir /code && cd /code \
    && git clone -b rel-2-2-7 https://gitlab.com/sherpa-team/sherpa.git && cd sherpa && autoreconf -fi \
    && ./configure --prefix=/usr/local \
      --enable-lhapdf=/usr/local --enable-fastjet=/usr/local \
      --enable-hepmc2=/usr/local --enable-openloops=/usr/local/OpenLoops \
      --enable-pyext --enable-gzip --enable-analysis \
      --enable-mpi CC=mpicc CXX=mpic++ \
    && make -j1 CXXFLAGS="-O2 -std=c++11" AM_CXXFLAGS="-O2 -std=c++11" \
    && make install \
    && mkdir /usr/local/Sherpa \
    && mv Examples /usr/local/Sherpa \
    && cd ../.. && rm -r /code
# wget https://sherpa.hepforge.org/downloads/?f=SHERPA-MC-2.2.9.tar.gz -O - | tar xf && cd SHERPA-MC-2.2.9

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
