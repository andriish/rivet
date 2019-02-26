FROM hepstore/rivet:2.7.0
MAINTAINER Andy Buckley <andy.buckley@cern.ch>

RUN dnf install -y boost-devel

RUN mkdir /code

# TODO: combine RUN rules to avoid the large intermediate layers

RUN cd /code && \
    wget https://www.hepforge.org/archive/lhapdf/LHAPDF-6.2.1.tar.gz -O- | tar xz && \
    cd LHAPDF-*/ && ./configure --prefix=/usr/local && \
    make -j5 && make install

#RUN lhapdf install MMHT2014{,n}lo68cl
RUN cd /usr/local/share/LHAPDF/ && \
    for pdf in MMHT2014{,n}lo68cl; do \
      wget https://lhapdf.hepforge.org/downloads?f=pdfsets/current/$pdf.tar.gz -O- | tar xz; \
    done

RUN cd /code && \
    wget https://www.hepforge.org/archive/thepeg/ThePEG-2.1.4.tar.bz2 -O- | tar xj && \
    cd ThePEG-*/ && ./configure --enable-shared --{prefix,with-{fastjet,hepmc,lhapdf,rivet}}=/usr/local && \
    make -j5 && make install

RUN cd /code && \
    wget https://www.hepforge.org/archive/herwig/Herwig-7.1.4.tar.bz2 -O- | tar xj && \
    cd Herwig-*/ && \
    ./configure --{prefix,with-{thepeg,fastjet}}=/usr/local && \
    make -j5 && make install

RUN rm -rf /code

WORKDIR /work
