FROM hepstore/rivet:3.1.2
LABEL maintainer="rivet-developers@cern.ch"

RUN mkdir /code && cd /code \
    && wget --no-verbose https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.3.0.tar.gz -O- | tar xz \
    && cd LHAPDF-*/ && ./configure --prefix=/usr/local \
    && make -j $(nproc --ignore=1) && make install \
    && cd ../.. && rm -rf /code

RUN lhapdf install MMHT2014{,n}lo68cl CT14{,n}lo

RUN mkdir /code && cd /code \
    && apt-get install -y mercurial libtool libboost-dev \
    && wget https://www.hepforge.org/archive/thepeg/ThePEG-2.2.1.tar.bz2 -O- | tar xj \
    && cd ThePEG-*/ && ./configure --enable-shared --{prefix,with-{fastjet,hepmc,lhapdf}}=/usr/local \
    && make -j $(nproc --ignore=1) && make install \
    && cd ../.. && rm -rf /code

RUN mkdir /code && cd /code \
    && wget https://www.hepforge.org/archive/herwig/Herwig-7.2.1.tar.bz2 -O- | tar xj \
    && cd Herwig-*/ \
    && ./configure --{prefix,with-{thepeg,fastjet}}=/usr/local \
    && make -j $(nproc --ignore=1) && make install \
    && cd ../.. && rm -rf /code

WORKDIR /work
