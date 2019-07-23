FROM hepstore/rivet:3.0.1
MAINTAINER Andy Buckley <andy.buckley@cern.ch>

RUN dnf install -y boost-devel

RUN mkdir /code && cd /code \
    && wget https://www.hepforge.org/archive/lhapdf/LHAPDF-6.2.3.tar.gz -O- | tar xz \
    && cd LHAPDF-*/ && ./configure --prefix=/usr/local \
    && make -j5 && make install \
    && cd ../.. && rm -rf /code

RUN lhapdf install MMHT2014{,n}lo68cl
# RUN cd /usr/local/share/LHAPDF/ && \
#     for pdf in MMHT2014{,n}lo68cl; do \
#       wget https://lhapdf.hepforge.org/downloads?f=pdfsets/current/$pdf.tar.gz -O- | tar xz; \
#     done

RUN mkdir /code && cd /code \
    && wget https://www.hepforge.org/archive/thepeg/ThePEG-2.1.5.tar.bz2 -O- | tar xj \
    && cd ThePEG-*/ && ./configure --enable-shared --{prefix,with-{fastjet,hepmc,lhapdf,rivet}}=/usr/local \
    && make -j5 && make install \
    && cd ../.. && rm -rf /code

RUN mkdir /code && cd /code \
    && wget https://www.hepforge.org/archive/herwig/Herwig-7.1.5.tar.bz2 -O- | tar xj \
    && cd Herwig-*/ \
    && ./configure --{prefix,with-{thepeg,fastjet}}=/usr/local \
    && make -j5 && make install \
    && cd ../.. && rm -rf /code

WORKDIR /work
