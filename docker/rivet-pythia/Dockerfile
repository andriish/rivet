FROM hepstore/rivet:3.1.2
LABEL maintainer="rivet-developers@cern.ch"

RUN mkdir /code && cd /code \
    && wget --no-verbose https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.3.0.tar.gz -O- | tar xz \
    && cd LHAPDF-*/ && ./configure --prefix=/usr/local \
    && make -j $(nproc --ignore=1) && make install \
    && cd ../.. && rm -r /code

RUN mkdir /code && cd /code \
    && apt-get -y install rsync && apt-get autoclean \
    && wget --no-verbose http://home.thep.lu.se/~torbjorn/pythia8/pythia8244.tgz -O- | tar xz \
    && cd pythia*/ && ./configure --enable-shared --{prefix,with-{hepmc2,lhapdf6}}=/usr/local \
    && make -j $(nproc --ignore=1) && make install \
    && cd ../.. && rm -r /code

RUN mkdir /code && cd /code \
    && wget --no-verbose https://agile.hepforge.org/downloads/?f=Sacrifice-1.1.2.tar.gz -O- | tar xz \
    && cd Sacrifice-*/ \
    && ./configure --prefix=/usr/local --with-pythia=/usr/local --with-hepmc=/usr/local --with-lhapdf=/usr/local \
    && make -j $(nproc --ignore=1) && make install \
    && cd ../.. && rm -r /code

WORKDIR /work
