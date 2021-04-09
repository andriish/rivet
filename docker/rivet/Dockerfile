ARG ARCH=ubuntu-gcc-hepmc3-py3
FROM hepstore/hepbase-${ARCH}-latex
LABEL maintainer="rivet-developers@cern.ch"
SHELL ["/bin/bash", "-c"]

ARG RIVET_VERSION
RUN mkdir /code && cd /code \
    && texhash \
    && wget --no-verbose https://rivet.hepforge.org/downloads/?f=Rivet-${RIVET_VERSION}.tar.gz -O Rivet-${RIVET_VERSION}.tar.gz \
    && tar xf Rivet-${RIVET_VERSION}.tar.gz && cd Rivet-${RIVET_VERSION}/ \
    && if [[ -d /usr/local/include/HepMC3 ]]; then HEPMC_FLAG="--with-hepmc3=/usr/local"; fi \
    && ./configure $HEPMC_FLAG && make -j $(nproc --ignore=1) && make install \
    && texhash \
    && echo "source /usr/local/etc/bash_completion.d/rivet-completion" > /etc/profile.d/rivet-completion.sh \
    && echo "source /usr/local/etc/bash_completion.d/yoda-completion" > /etc/profile.d/yoda-completion.sh \
    && texhash \
    && rm -rf /code

ADD bash.bashrc /etc/

WORKDIR /work
