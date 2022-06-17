ARG ARCH=ubuntu-gcc-hepmc3-py3
FROM hepstore/hepbase-${ARCH}-latex
LABEL maintainer="rivet-developers@cern.ch"
SHELL ["/bin/bash", "--login", "-c"]

ARG YODA_BRANCH
ARG RIVET_BRANCH
RUN mkdir /code && cd /code \
    && texhash \
    && git clone --depth 1 --branch $YODA_BRANCH https://gitlab.com/hepcedar/yoda.git yoda/ \
    && cd yoda/ \
    && autoreconf -i \
    && ./configure && make -j $(nproc --ignore=1) && make install \
    && cd .. \
    && git clone --depth 1 --branch $RIVET_BRANCH https://gitlab.com/hepcedar/rivet.git rivet/ \
    && cd rivet/ \
    && autoreconf -i \
    && if [[ -d /usr/local/include/HepMC3 ]]; then HEPMC_FLAG="--with-hepmc3=/usr/local"; fi \
    && ./configure $HEPMC_FLAG && make -j $(nproc --ignore=1) && make install \
    && texhash \
    && echo "source /usr/local/etc/bash_completion.d/rivet-completion" > /etc/profile.d/rivet-completion.sh \
    && echo "source /usr/local/etc/bash_completion.d/yoda-completion" > /etc/profile.d/yoda-completion.sh \
    && texhash \
    && rm -rf /code

# TODO: remove when hepbase is fixed
RUN export DEBIAN_FRONTEND=noninteractive \
    && PYV=$(python -c 'from __future__ import print_function; import sys; print("{i.major}.{i.minor}".format(i=sys.version_info))') \
    && echo "PYTHONPATH=\"\$PYTHONPATH:/usr/local/lib/python$PYV/site-packages/\"" >> /etc/profile.d/05-usrlocal.sh

#ADD bash.bashrc /etc/

WORKDIR /work
