ARG BASE=ubuntu-gcc-hepmc3-py3-latex
FROM hepstore/hepbase-${BASE}
LABEL maintainer="rivet-developers@cern.ch"

ARG RIVET_VERSION=3.1.2
RUN mkdir /code && cd /code \
    && wget --no-verbose https://rivet.hepforge.org/downloads/?f=Rivet-${RIVET_VERSION}.tar.gz -O Rivet-${RIVET_VERSION}.tar.gz \
    && tar xf Rivet-${RIVET_VERSION}.tar.gz && cd Rivet-${RIVET_VERSION}/ \
    && test -d /usr/local/include/HepMC3 && HEPMC_FLAG="--with-hepmc3=/usr/local" \
    && ./configure $HEPMC_FLAG && make -j $(nproc --ignore=1) && make install \
    && echo "source /usr/local/share/Rivet/rivet-completion" > /etc/profile.d/rivet-completion.sh \
    && echo "source /usr/local/share/YODA/yoda-completion" > /etc/profile.d/yoda-completion.sh \
    && texhash \
    && rm -rf /code

WORKDIR /work
