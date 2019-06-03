FROM hepstore/rivet:2.7.2
MAINTAINER Andy Buckley <andy.buckley@cern.ch>

RUN dnf install -y rsync

RUN mkdir /code && cd /code \
    && wget https://www.hepforge.org/archive/lhapdf/LHAPDF-6.2.3.tar.gz -O- | tar xz \
    && cd LHAPDF-*/ && ./configure --prefix=/usr/local \
    && make -j5 && make install \
    && cd ../.. && rm -r /code

RUN mkdir /code && cd /code \
    && wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8240.tgz -O- | tar xz \
    && cd pythia*/ && ./configure --enable-shared --{prefix,with-{hepmc2,lhapdf6}}=/usr/local \
    && make -j5 && make install \
    && cd ../.. && rm -r /code

RUN mkdir /code && cd /code \
    && wget https://agile.hepforge.org/downloads/?f=Sacrifice-1.1.2.tar.gz -O- | tar xz \
    && cd Sacrifice-*/ \
    && ./configure --{prefix,with-{pythia,hepmc,LHAPDF}}=/usr/local \
    && make -j5 && make install \
    && cd ../.. && rm -r /code

WORKDIR /work
