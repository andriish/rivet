FROM hepstore/hepbase-ubuntu-gcc-hepmc2-py3-latex
LABEL maintainer="rivet-developers@cern.ch"

ARG RIVET_VERSION=3.1.1

RUN mkdir /code && cd /code \
    && wget https://rivet.hepforge.org/downloads/?f=Rivet-${RIVET_VERSION}.tar.gz -O Rivet-${RIVET_VERSION}.tar.gz \
    && tar xf Rivet-${RIVET_VERSION}.tar.gz && cd Rivet-${RIVET_VERSION}/ \
    && ./configure && make -j4 && make install \
    && echo "source /usr/local/share/Rivet/rivet-completion" > /etc/profile.d/rivet-completion.sh \
    && echo "source /usr/local/share/YODA/yoda-completion" > /etc/profile.d/yoda-completion.sh \
    && texhash \
    && rm -rf /code

WORKDIR /work
